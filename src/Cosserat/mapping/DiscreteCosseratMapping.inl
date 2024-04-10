/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once
#include "DiscreteCosseratMapping.h"

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/gl/template.h>
#include <sofa/helper/visual/DrawTool.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/Quat.h>

#include <string>


namespace sofa::component::mapping
{

    using sofa::core::objectmodel::BaseContext ;
    using sofa::helper::AdvancedTimer;
    using sofa::helper::WriteAccessor;
    using sofa::defaulttype::SolidTypes ;
    using sofa::type::RGBAColor;


template <class TIn1, class TIn2, class TOut>
DiscreteCosseratMapping<TIn1, TIn2, TOut>::DiscreteCosseratMapping()
    : m_fromModel1(NULL)
    , m_fromModel2(NULL)
    , m_toModel(NULL)
    , d_deformationAxis(initData(&d_deformationAxis, (int)1, "deformationAxis",
                                 "the axis in which we want to show the deformation.\n"))
    , d_max(initData(&d_max, (Real2)1.0e-2, "max",
                                 "the maximum of the deformation.\n"))
    , d_min(initData(&d_min, (Real2)0.0, "min",
                                 "the minimum of the deformation.\n"))
    , d_radius(initData(&d_radius, (Real2)0.05, "radius",
                                 "the axis in which we want to show the deformation.\n"))
    , d_drawMapBeam(initData(&d_drawMapBeam, true,"nonColored", "if this parameter is false, you draw the beam with "
                                                                "color according to the force apply to each beam"))
    , d_color(initData(&d_color, sofa::type::RGBAColor(40/255.0, 104/255.0, 137/255.0, 0.8) ,"color", "The default beam color"))
    , d_index(initData(&d_index, "index", "if this parameter is false, you draw the beam with color "
                                                          "according to the force apply to each beam"))
    , d_baseIndex(initData(&d_baseIndex, (unsigned int) 0, "baseIndex", "This parameter defines the index of the rigid "
                                                                        "base of Cosserat models, 0 by default this can"
                                                                        "take another value if the rigid base is given "
                                                                        "by another body."))
{}

template <class TIn1, class TIn2, class TOut>
void DiscreteCosseratMapping<TIn1, TIn2, TOut>::init()
{
    if(this->getFromModels1().empty() || this->getFromModels2().empty() || this->getToModels().empty())
    {
        msg_error() << "Error while initializing ; input getFromModels1/getFromModels2/output not found" ;
        return;
    }

    reinit();
    Inherit1::init();

    m_colorMap.setColorScheme("Blue to Red");
    m_colorMap.reinit();
}

template <class TIn1, class TIn2, class TOut>
void DiscreteCosseratMapping<TIn1, TIn2, TOut>::reinit()
{
  m_fromModel1 = this->getFromModels1()[0]; // Cosserat deformations (torsion and bending), in local frame
  m_fromModel2 = this->getFromModels2()[0]; // Cosserat base, in global frame
  m_toModel = this->getToModels()[0];  // Cosserat rigid frames, in global frame

  // Fill the initial vectorvoid init() override;
  const OutDataVecCoord* xFromData = m_toModel->read(core::ConstVecCoordId::position());
  const OutVecCoord xFrom = xFromData->getValue();

  m_vecTransform.clear();
  for (unsigned int i = 0; i < xFrom.size(); i++) {
    m_vecTransform.push_back(xFrom[i]);
  }

  if(d_debug.getValue())
    msg_info("DiscreteCosseratMapping")<< " m_vecTransform : "<< m_vecTransform;

  this->initialize();
}

template <class TIn1, class TIn2, class TOut>
void DiscreteCosseratMapping<TIn1, TIn2, TOut>::apply(
        const core::MechanicalParams* /* mparams */, const type::vector<OutDataVecCoord*>& dataVecOutPos,
        const type::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const type::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{
    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    ///Do Apply
    //We need only one input In model and input Root model (if present)
    const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();

    const auto sz = d_curv_abs_frames.getValue().size();
    OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
    out.resize(sz);
    const auto baseIndex = d_baseIndex.getValue();

    // update the Exponential matrices according to new deformation
    // Here we update m_framesExponentialSE3Vectors & m_nodesExponentialSE3Vectors
    this->update_ExponentialSE3(in1);
    /* from cossserat to SOFA frame*/
    Transform frame0 = Transform(In2::getCPos(in2[baseIndex]),In2::getCRot(in2[baseIndex]));
    for(unsigned int i=0; i<sz; i++){
        Transform frame = frame0;
        for (unsigned int u = 0; u < m_indicesVectors[i]; u++) {
            frame *= m_nodesExponentialSE3Vectors[u];
        }
        frame *= m_framesExponentialSE3Vectors[i];

        type::Vec3 v = frame.getOrigin();
        type::Quat q = frame.getOrientation();
        out[i] = OutCoord(v,q);
    }
    // @todo
    m_index_input = 0;
    dataVecOutPos[0]->endEdit();
}


//template <class TIn1, class TIn2, class TOut>
//Matrix4 DiscreteCosseratMapping<TIn1, TIn2, TOut>::computeLogarithm(const double & x, const defaulttype::Matrix4 &gX){

//    // Compute theta before everything
//    const double theta = computeTheta(x, gX);
//    Matrix4 I4; I4.clear(); I4.identity();
//    Matrix4 log_gX; log_gX.clear();


//    double csc_theta = 1.0/(sin(x * theta/2.0));
//    double sec_theta = 1.0/(cos(x * theta/2.0));
//    double cst = (1.0/8) * (csc_theta*csc_theta*csc_theta) * sec_theta;
//    double x_theta = x*theta;
//    double cos_2Xtheta = cos(2.0 * x_theta);
//    double cos_Xtheta = cos(x_theta);
//    double sin_2Xtheta = sin(2.0 *x_theta);
//    double sin_Xtheta = sin(x_theta);

//    if(theta <= std::numeric_limits<double>::epsilon()) log_gX = I4;
//    else {
//        log_gX  = cst * ((x_theta*cos_2Xtheta - sin_Xtheta)*I4 -
//                         (x_theta*cos_Xtheta + 2.0*x_theta*cos_2Xtheta - sin_Xtheta -sin_2Xtheta)*gX +
//                         (2.0*x_theta*cos_Xtheta + x_theta*cos_2Xtheta-sin_Xtheta - sin_2Xtheta) *(gX*gX)-
//                         (x_theta*cos_Xtheta - sin_Xtheta)*(gX*gX*gX));
//    }

//    return log_gX;
//}

//template<class In1VecCoord, class Mat6x6>
//void computeViolation(In1VecCoord& inDeform, const helper::vector<double> m_framesLengthVectors, const
//                       size_t sz, std::function<double(int i, int j)> f)
//{


//    for (std::size_t i = 0; i < sz; i++) {
//        Mat6x6 temp ;

//        type::Vec3 k = inDeform[m_indicesVectors[i]-1];
//        double  x = m_framesLengthVectors[i];
//        compute_Tang_Exp(x,k,temp) ;
//        m_framesTangExpVectors.push_back(temp);

////        if (d_debug.getValue()){
////            printf("__________________________________________\n");
////            std::cout << "x :"<< x << "; k :"<< k << std::endl;
////            std::cout<< "m_framesTangExpVectors :"<< m_framesTangExpVectors[i] << std::endl;
////        }
//    }
//}


template <class TIn1, class TIn2, class TOut>
void DiscreteCosseratMapping<TIn1, TIn2, TOut>:: applyJ(
        const core::MechanicalParams* /* mparams */, const type::vector< OutDataVecDeriv*>& dataVecOutVel,
        const type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const type::vector<const In2DataVecDeriv*>& dataVecIn2Vel) {

    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
        return;
    const In1VecDeriv& in1 = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& in2_vecDeriv = dataVecIn2Vel[0]->getValue();
    OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();
    const auto baseIndex = d_baseIndex.getValue();

    // This is the vector of the curv abscissa, X in the paper
    helper::ReadAccessor<Data<type::vector<double>>> curv_abs_section =  d_curv_abs_section;
    helper::ReadAccessor<Data<type::vector<double>>> curv_abs_frames = d_curv_abs_frames;

    // Compute the tangent Exponential SE3 vectors
    const In1VecCoord& inDeform = m_fromModel1->read(core::ConstVecCoordId::position())->getValue();
    this->update_TangExpSE3(inDeform, curv_abs_section.ref(), curv_abs_frames.ref());

    //Get base velocity as input this is also called eta
    m_nodesVelocityVectors.clear();
    Deriv2 _baseVelocity;
    if (!in2_vecDeriv.empty())
        _baseVelocity = in2_vecDeriv[baseIndex];
    //convert to Vec6
    type::Vec6 baseVelocity;
    for (auto u=0; u<6; u++) {baseVelocity[u] = _baseVelocity[u];}

    //Apply the local transform i.e from SOFA's frame to Cosserat's frame
    const In2VecCoord& xfrom2Data = m_fromModel2->read(core::ConstVecCoordId::position())->getValue();
    Transform TInverse = Transform(xfrom2Data[baseIndex].getCenter(), xfrom2Data[baseIndex].getOrientation()).inversed();
    Mat6x6 P = this->build_projector(TInverse);
    type::Vec6 baseLocalVelocity = P * baseVelocity;
    m_nodesVelocityVectors.push_back(baseLocalVelocity);
    if(d_debug.getValue())
        std::cout << "Base local Velocity :"<< baseLocalVelocity <<std::endl;

    //Compute velocity at nodes
    for (unsigned int i = 1 ; i < curv_abs_section.size(); i++) {
        Transform Trans = m_nodesExponentialSE3Vectors[i].inversed();
        Mat6x6 Adjoint; Adjoint.clear();
        this->computeAdjoint(Trans, Adjoint);

        type::Vec6 Xi_dot = Vec6(in1[i-1],type::Vec3(0.0,0.0,0.0)) ;
        Vec6 temp = Adjoint * (m_nodesVelocityVectors[i-1] + m_nodesTangExpVectors[i] * Xi_dot );
        m_nodesVelocityVectors.push_back(temp);
        if(d_debug.getValue())
            std::cout<< "Node velocity : "<< i << " = " << temp<< std::endl;
    }

    const OutVecCoord& out = m_toModel->read(core::ConstVecCoordId::position())->getValue();
    auto sz = curv_abs_frames.size();
    outVel.resize(sz);
    for (unsigned int i = 0 ; i < sz; i++) {
        Transform Trans = m_framesExponentialSE3Vectors[i].inversed();
        Mat6x6 Adjoint; Adjoint.clear();
        this->computeAdjoint(Trans, Adjoint);

        type::Vec6 Xi_dot = Vec6(in1[m_indicesVectors[i]-1],type::Vec3(0.0,0.0,0.0)) ;
        Vec6 temp = Adjoint * (m_nodesVelocityVectors[m_indicesVectors[i]-1] + m_framesTangExpVectors[i] * Xi_dot ); // eta

        auto T = Transform(out[i].getCenter(), out[i].getOrientation());
        Mat6x6 Proj = this->build_projector(T);

        outVel[i] = Proj * temp;

        if(d_debug.getValue())
            std::cout<< "Frame velocity : "<< i << " = " << temp<< std::endl;
    }
    dataVecOutVel[0]->endEdit();
    m_index_input = 0;
}


template <class TIn1, class TIn2, class TOut>
void DiscreteCosseratMapping<TIn1, TIn2, TOut>:: applyJT(
        const core::MechanicalParams* /*mparams*/, const type::vector< In1DataVecDeriv*>& dataVecOut1Force,
        const type::vector< In2DataVecDeriv*>& dataVecOut2Force,
        const type::vector<const OutDataVecDeriv*>& dataVecInForce)  {

    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

    const OutVecDeriv& in = dataVecInForce[0]->getValue();

    In1VecDeriv& out1 = *dataVecOut1Force[0]->beginEdit();
    In2VecDeriv& out2 = *dataVecOut2Force[0]->beginEdit();
    const auto baseIndex = d_baseIndex.getValue();

    const OutVecCoord& frame = m_toModel->read(core::ConstVecCoordId::position())->getValue();
    const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
    const In1VecCoord x1from = x1fromData->getValue();
    type::vector<Vec6> local_F_Vec;  local_F_Vec.clear();

    out1.resize(x1from.size());

    //convert the input from Deriv type to vec6 type, for the purpose of the matrix vector multiplication
    for (unsigned int var = 0; var < in.size(); ++var) {
        type::Vec6 vec;
        for(unsigned j = 0; j < 6; j++) vec[j] = in[var][j];
        //Convert input from global frame(SOFA) to local frame
        Transform _T = Transform(frame[var].getCenter(),frame[var].getOrientation());
        Mat6x6 P_trans =(this->build_projector(_T)); P_trans.transpose();
        type::Vec6 local_F = P_trans * vec;
        local_F_Vec.push_back(local_F);
    }

    //Compute output forces
    auto sz = m_indicesVectors.size();
    auto index =  m_indicesVectors[sz-1];
    m_totalBeamForceVectors.clear();
    m_totalBeamForceVectors.resize(sz);

    Vec6 F_tot; F_tot.clear();
    m_totalBeamForceVectors.push_back(F_tot);

    Mat3x6 matB_trans; matB_trans.clear();
    for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;


    for (auto s = sz ; s-- ; ) {
        Mat6x6 coAdjoint;

        this->compute_coAdjoint(m_framesExponentialSE3Vectors[s], coAdjoint);  // m_framesExponentialSE3Vectors[s] computed in apply
        Vec6 node_F_Vec = coAdjoint * local_F_Vec[s];
        Mat6x6 temp = m_framesTangExpVectors[s];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
        temp.transpose();
        type::Vec3 f = matB_trans * temp * node_F_Vec;

        if(index != m_indicesVectors[s]){
            index--;
            //bring F_tot to the reference of the new beam
            this->compute_coAdjoint(m_nodesExponentialSE3Vectors[index],coAdjoint);  //m_nodesExponentialSE3Vectors computed in apply
            F_tot = coAdjoint * F_tot;
            Mat6x6 temp = m_nodesTangExpVectors[index];
            temp.transpose();
            //apply F_tot to the new beam
            type::Vec3 temp_f = matB_trans * temp * F_tot;
            out1[index-1] += temp_f;
        }
        if(d_debug.getValue())
            std::cout << "f at s ="<< s <<" and index"<< index <<  " is : "<< f << std::endl;


        //compute F_tot
        F_tot += node_F_Vec;
        out1[m_indicesVectors[s]-1] += f;
    }

    Transform frame0 = Transform(frame[0].getCenter(),frame[0].getOrientation());
    Mat6x6 M = this->build_projector(frame0);
    out2[baseIndex] += M * F_tot;

    if(d_debug.getValue()){
        std::cout << "Node forces "<< out1 << std::endl;
        std::cout << "base Force: "<< out2[baseIndex] << std::endl;
    }

    dataVecOut1Force[0]->endEdit();
    dataVecOut2Force[0]->endEdit();
}

//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void DiscreteCosseratMapping<TIn1, TIn2, TOut>::applyJT(
        const core::ConstraintParams*/*cparams*/ , const type::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
        const type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
        const type::vector<const OutDataMatrixDeriv*>& dataMatInConst)
{
    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
        return;

    //We need only one input In model and input Root model (if present)
    In1MatrixDeriv& out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the strain space (reduced coordinate)
    In2MatrixDeriv& out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the reference frame (base frame)
    const OutMatrixDeriv& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames

    const OutVecCoord& frame = m_toModel->read(core::ConstVecCoordId::position())->getValue();
    const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
    const In1VecCoord x1from = x1fromData->getValue();

    Mat3x6 matB_trans; matB_trans.clear();
    for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;

    type::vector< std::tuple<int,Vec6> > NodesInvolved;
    type::vector< std::tuple<int,Vec6> > NodesInvolvedCompressed;
    //helper::vector<Vec6> NodesConstraintDirection;

    typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        if (d_debug.getValue()){
            std::cout<<"************* Apply JT (MatrixDeriv) iteration on line ";
            std::cout<<rowIt.index();
            std::cout<<"*************  "<<std::endl;
        }
        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
        typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

        // Creates a constraints if the input constraint is not empty.
        if (colIt == colItEnd)
        {
            if (d_debug.getValue()){
                std::cout<<"no column for this constraint"<<std::endl;
            }
            continue;
        }
        typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
        typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());

        NodesInvolved.clear();
        while (colIt != colItEnd)
        {
            int childIndex = colIt.index();

            const OutDeriv valueConst_ = colIt.val();
            type::Vec6 valueConst;
            for(unsigned j = 0; j < 6; j++) valueConst[j] = valueConst_[j];

            int indexBeam =  m_indicesVectors[childIndex];

            Transform _T = Transform(frame[childIndex].getCenter(),frame[childIndex].getOrientation());
            Mat6x6 P_trans =(this->build_projector(_T));
            P_trans.transpose();

            Mat6x6 coAdjoint;
            this->compute_coAdjoint(m_framesExponentialSE3Vectors[childIndex], coAdjoint);  // m_framesExponentialSE3Vectors[s] computed in apply
            Mat6x6 temp = m_framesTangExpVectors[childIndex];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
            temp.transpose();

            type::Vec6 local_F =  coAdjoint * P_trans * valueConst; // constraint direction in local frame of the beam.

            type::Vec3 f = matB_trans * temp * local_F; // constraint direction in the strain space.

            o1.addCol(indexBeam-1, f);
            std::tuple<int,Vec6> test = std::make_tuple(indexBeam, local_F);

            NodesInvolved.push_back(test);
            colIt++;

        }
        if (d_debug.getValue()){
            std::cout<<"==> NodesInvolved : "<<std::endl;
            for (size_t i = 0; i < NodesInvolved.size(); i++)
                std::cout << "index :" <<get<0>(NodesInvolved[i]) << " force :"
                          << get<1>(NodesInvolved[i]) << "\n ";
        }

        // sort the Nodes Invoved by decreasing order
        std::sort(begin(NodesInvolved), end(NodesInvolved),
                  [](std::tuple<int, Vec6> const &t1, std::tuple<int, Vec6> const &t2) {
            return std::get<0>(t1) > std::get<0>(t2); // custom compare function
        } );

        NodesInvolvedCompressed.clear();

        for (unsigned n=0; n<NodesInvolved.size(); n++)
        {
            std::tuple<int,Vec6> test_i = NodesInvolved[n];
            int numNode_i= std::get<0>(test_i);
            Vec6 cumulativeF =std::get<1>(test_i);

            if (n<NodesInvolved.size()-1)
            {
                std::tuple<int,Vec6> test_i1 = NodesInvolved[n+1];
                int numNode_i1= std::get<0>(test_i1);

                while (numNode_i == numNode_i1)
                {
                    cumulativeF += std::get<1>(test_i1);
                    //// This was if ((n!=NodesInvolved.size()-2)||(n==0)) before and I change it to
                    /// if ((n!=NodesInvolved.size()-1)||(n==0)) since the code can't leave the will loop
                    if ((n!=NodesInvolved.size()-1)||(n==0)){
                        n++;
                        break;
                    }
                    test_i1 = NodesInvolved[n+1];
                    numNode_i1= std::get<0>(test_i1);
                }

            }
            NodesInvolvedCompressed.push_back(std::make_tuple(numNode_i, cumulativeF));
        }

        if (d_debug.getValue()){
            std::cout<<" NodesInvolved after sort and compress"<<std::endl;
            for (size_t i = 0; i < NodesInvolvedCompressed.size(); i++)
                std::cout << "index :" <<get<0>(NodesInvolvedCompressed[i]) << " force :"
                          << get<1>(NodesInvolvedCompressed[i]) << "\n ";
        }

        for (unsigned n=0; n<NodesInvolvedCompressed.size(); n++)
        {
            std::tuple<int,Vec6> test = NodesInvolvedCompressed[n];
            int numNode= std::get<0>(test);
            int i = numNode;
            Vec6 CumulativeF = std::get<1>(test);

            while(i>0)
            {
                //cumulate on beam frame
                Mat6x6 coAdjoint;
                this->compute_coAdjoint(m_nodesExponentialSE3Vectors[i-1],coAdjoint);  //m_nodesExponentialSE3Vectors computed in apply
                CumulativeF = coAdjoint * CumulativeF;
                // transfer to strain space (local coordinates)
                Mat6x6 temp = m_nodesTangExpVectors[i-1];
                temp.transpose();
                type::Vec3 temp_f = matB_trans * temp * CumulativeF;

                if(i>1) o1.addCol(i-2, temp_f);
                i--;
            }

            Transform frame0 = Transform(frame[0].getCenter(),frame[0].getOrientation());
            Mat6x6 M = this->build_projector(frame0);

            Vec6 base_force = M * CumulativeF;
            o2.addCol(d_baseIndex.getValue(), base_force);
        }
    }

    //"""END ARTICULATION SYSTEM MAPPING"""
    dataMatOut1Const[0]->endEdit();
    dataMatOut2Const[0]->endEdit();
}

template <class TIn1, class TIn2, class TOut>
void DiscreteCosseratMapping<TIn1, TIn2, TOut>::computeBBox(const core::ExecParams*, bool)
{
//  const VecCoord& x = getVertices(); //m_vertices.getValue();
  const OutVecCoord& x = m_toModel->read(core::ConstVecCoordId::position())->getValue();

  SReal minBBox[3] = {std::numeric_limits<Real1>::max(),std::numeric_limits<Real1>::max(),std::numeric_limits<Real1>::max()};
  SReal maxBBox[3] = {-std::numeric_limits<Real1>::max(),-std::numeric_limits<Real1>::max(),-std::numeric_limits<Real1>::max()};
  for (std::size_t i = 0; i < x.size(); i++)
  {
    const OutCoord& p = x[i];
    for (int c=0; c<3; c++)
    {
      if (p[c] > maxBBox[c]) maxBBox[c] = p[c];
      if (p[c] < minBBox[c]) minBBox[c] = p[c];
    }
  }
  this->f_bbox.setValue(sofa::type::TBoundingBox<SReal>(minBBox,maxBBox));
}


template <class TIn1, class TIn2, class TOut>
void DiscreteCosseratMapping<TIn1, TIn2, TOut>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMechanicalMappings()) return;

    // draw cable
    typedef RGBAColor RGBAColor;

    const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();

    const OutDataVecCoord* xfromData = m_toModel->read(core::ConstVecCoordId::position());
    const OutVecCoord xData = xfromData->getValue();
    type::vector<type::Vec3> positions;
    type::vector<sofa::type::Quat<Real2>> Orientation;
    positions.clear();
    Orientation.clear();
    unsigned int sz = xData.size();
    for (unsigned int i = 0; i<sz; i++){
        positions.push_back(xData[i].getCenter());
        Orientation.push_back(xData[i].getOrientation());
    }

    //Get access articulated
    const In1DataVecCoord* artiData = m_fromModel1->read(core::ConstVecCoordId::position());
    const In1VecCoord xPos = artiData->getValue();

    RGBAColor drawColor = d_color.getValue();
    for (unsigned int i=0; i<sz-1; i++)
        vparams->drawTool()->drawCylinder(positions[i], positions[i+1], d_radius.getValue(), drawColor);

    //Define color map
    Real2 min = d_min.getValue();
    Real2 max = d_max.getValue();
    helper::ColorMap::evaluator<Real2> _eval = m_colorMap.getEvaluator(min, max);

    glLineWidth(d_radius.getValue());
    glBegin(GL_LINES);
    if(d_drawMapBeam.getValue()){
        sofa::type::RGBAColor _color = d_color.getValue();
        RGBAColor colorL = RGBAColor(_color[0],_color[1],_color[2],_color[3]);
        glColor4f(colorL[0], colorL[1], colorL[2],colorL[3]);
        for (unsigned int i=0; i<sz-1; i++) {
            vparams->drawTool()->drawLine(positions[i],positions[i+1],colorL);
        }
    }else {
        int j = 0;
        type::vector<int> index = d_index.getValue();
        for (unsigned int i=0; i<sz-1; i++) {
            j = m_indicesVectorsDraw[i]-1; // to get the articulation on which the frame is related to
            RGBAColor color =  _eval(xPos[j][d_deformationAxis.getValue()]);
            vparams->drawTool()->drawLine(positions[i],positions[i+1],color);
        }
    }
    glLineWidth(1);
    if (!vparams->displayFlags().getShowMappings())
        if(!d_debug.getValue()) return;
    glEnd();
}

} // namespace sofa
