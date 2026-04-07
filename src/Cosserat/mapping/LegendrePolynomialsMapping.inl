//
// Created by younes on 17/11/2021.
//
#pragma once
#include <Cosserat/config.h>
#include <Cosserat/mapping/LegendrePolynomialsMapping.h>

#include <sofa/component/mapping/nonlinear/RigidMapping.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/helper/io/Mesh.h>
#include <sofa/helper/io/SphereLoader.h>

namespace Cosserat::mapping {

	template<class TIn, class TOut>
	LegendrePolynomialsMapping<TIn, TOut>::LegendrePolynomialsMapping() :
		Inherit(), 
        d_order(initData(&d_order, (unsigned) 3, "order", "The order of Legendre polynomials")),
		d_vectorOfCurvilinearAbscissa(initData(&d_vectorOfCurvilinearAbscissa, "curvAbscissa",
											   "Vector of curvilinear Abscissa element of [0, 1]")) {}

	template<class TIn, class TOut>
	double LegendrePolynomialsMapping<TIn, TOut>::legendrePoly(unsigned int n, const double x) {
		// Boost provides a highly optimized and mathematically correct implementation
		// of Legendre polynomials, preventing the O(2^n) recursion and Bonnet formula errors.
		return boost::math::legendre_p(n, x);
	}

	template<class TIn, class TOut>
	void LegendrePolynomialsMapping<TIn, TOut>::reinit() {
		m_matOfCoeffs.clear();
		auto curvAbs = d_vectorOfCurvilinearAbscissa.getValue();
		auto sz = curvAbs.size();

		// Note: loop deliberately starts at i=1, skipping curvAbs[0].
		// This is kept for backward compatibility as existing scenes might rely on this shift.
		for (unsigned int i = 1; i < sz; i++) {
			vector<double> coeffsOf_i;
			coeffsOf_i.clear();
			for (unsigned int order = 0; order < d_order.getValue(); order++)
				coeffsOf_i.push_back(legendrePoly(order, curvAbs[i]));

			m_matOfCoeffs.push_back(coeffsOf_i);
		}
	}


	template<class TIn, class TOut>
	void LegendrePolynomialsMapping<TIn, TOut>::init() {
		// Compute the coefficients for each curv_abs at all orders of the polynomials
		reinit();

		Inherit1::init();
	}


	template<class TIn, class TOut>
	void LegendrePolynomialsMapping<TIn, TOut>::apply(const sofa::core::MechanicalParams * /*mparams*/,
													  Data<VecCoord> &dOut, const Data<InVecCoord> &dIn) {
		sofa::helper::ReadAccessor<Data<InVecCoord>> in = dIn;
		sofa::helper::WriteOnlyAccessor<Data<VecCoord>> out = dOut;
		const auto sz = d_vectorOfCurvilinearAbscissa.getValue().size();
		out.resize(sz - 1);

		for (unsigned int i = 0; i < sz - 1; i++) {
			Vec3 Xi(0.0, 0.0, 0.0);
			for (unsigned int j = 0; j < in.size(); j++)
				Xi += m_matOfCoeffs[i][j] * in[j];

			out[i] = Xi;
		}
	}

	template<class TIn, class TOut>
	void LegendrePolynomialsMapping<TIn, TOut>::applyJ(const sofa::core::MechanicalParams * /*mparams*/,
													   Data<VecDeriv> &dOut, const Data<InVecDeriv> &dIn) {
		sofa::helper::WriteOnlyAccessor<Data<VecDeriv>> velOut = dOut;
		sofa::helper::ReadAccessor<Data<InVecDeriv>> velIn = dIn;

		const auto sz = d_vectorOfCurvilinearAbscissa.getValue().size();
		velOut.resize(sz - 1);
		for (sofa::Index i = 0; i < sz - 1; ++i) {
			Vec3 vel(0.0, 0.0, 0.0);
			for (unsigned int j = 0; j < velIn.size(); j++)
				vel += m_matOfCoeffs[i][j] * velIn[j];
			velOut[i] = vel;
		}
	}

	template<class TIn, class TOut>
	void LegendrePolynomialsMapping<TIn, TOut>::applyJT(const sofa::core::MechanicalParams * /*mparams*/,
														Data<InVecDeriv> &dOut, const Data<VecDeriv> &dIn) {
		sofa::helper::WriteAccessor<Data<InVecDeriv>> out = dOut;
		sofa::helper::ReadAccessor<Data<VecDeriv>> in = dIn;
		const unsigned int numDofs = this->getFromModel()->getSize();
		out.resize(numDofs);
		for (unsigned int cI = 0; cI < out.size(); cI++) {
			// Clamp coefficient index to avoid out-of-bounds access if numDofs != d_order
			const unsigned int max_col = d_order.getValue();
			if (cI >= max_col)
				continue;

			for (sofa::Index i = 0; i < in.size(); ++i) {
				// std::cout << " cI:" << cI << " i:"<< i <<" m_matOfCoeffs[i][cI] : "<< m_matOfCoeffs[i][cI] * in[i]<<
				// std::endl;
				//@todo use alpha factor
				out[cI] += m_matOfCoeffs[i][cI] * in[i];
			}
		}
	}

	// Propagate the constraint through the Legendre polynomials mapping
	template<class TIn, class TOut>
	void LegendrePolynomialsMapping<TIn, TOut>::applyJT(const sofa::core::ConstraintParams * /*cparams*/,
														Data<InMatrixDeriv> &dOut, const Data<OutMatrixDeriv> &dIn) {
		InMatrixDeriv &out = *dOut.beginEdit();
		const OutMatrixDeriv &in = dIn.getValue();

		const unsigned int numDofs = this->getFromModel()->getSize();
		const unsigned int max_col = d_order.getValue();

		typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

		for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt) {
			typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
			typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

			if (colIt == colItEnd)
				continue;

			typename InMatrixDeriv::RowIterator o = out.writeLine(rowIt.index()); // we store the constraint number
			while (colIt != colItEnd) {
				int childIndex = colIt.index();
				const OutDeriv f_It = colIt.val();
				for (unsigned int order = 0; order < numDofs; order++) {
					if (order >= max_col)
						continue;

					InDeriv f;
					f = m_matOfCoeffs[childIndex][order] * f_It;
					o.addCol(order, f);
				}
				colIt++;
			}
		}

		dOut.endEdit();
	}

} // namespace Cosserat::mapping
