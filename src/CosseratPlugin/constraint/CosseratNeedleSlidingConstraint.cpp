
/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* This component is not open-source                                           *
*                                                                             *
* Authors: Yinoussa Adagolodjo/adagolodjo@protonamil.com                      *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define SOFA_COMPONENT_CONSTRAINTSET_COSSERATNEEDLESLIDINGCONSTRAINT_CPP

#include<sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>
#include "CosseratNeedleSlidingConstraint.inl"


namespace sofa::component::constraintset
{

using sofa::defaulttype::Rigid3Types;
using namespace sofa::helper;
using namespace sofa::core;

////////////////////////////////////////////    FACTORY    //////////////////////////////////////////////
// Registering the component
// see: http://wiki.sofa-framework.org/wiki/ObjectFactory
// 1-RegisterObject("description") + .add<> : Register the component
// 2-.add<>(true) : Set default template

class RigidImpl {};


template<>
class CosseratNeedleSlidingConstraintSpecialization<RigidImpl>
{
public:

    template <class T>
    static void buildConstraintMatrix(CosseratNeedleSlidingConstraint<T>& self,
                                      const ConstraintParams* cParams,
                                      typename CosseratNeedleSlidingConstraint<T>::DataMatrixDeriv &cMatrix,
                                      unsigned int &cIndex,
                                      const typename CosseratNeedleSlidingConstraint<T>::DataVecCoord &x)
    {
        if(self.d_componentState.getValue() != ComponentState::Valid)
            return ;

        SOFA_UNUSED(cParams);

        typename CosseratNeedleSlidingConstraint<T>::MatrixDeriv& matrix = *cMatrix.beginEdit();
        typename CosseratNeedleSlidingConstraint<T>::VecCoord positions = x.getValue();

        const type::Vec<3, typename CosseratNeedleSlidingConstraint<T>::Real> cy(0,1,0), cz(0,0,1);
        const type::Vec<3, typename CosseratNeedleSlidingConstraint<T>::Real> vZero(0,0,0);


        self.m_constraintId= cIndex;

        type::Vec<6,bool> use = self.d_useDirections.getValue();

        for (unsigned int i=0; i<positions.size(); i++)
        {
            if (use[1]){
              typename CosseratNeedleSlidingConstraint<T>::MatrixDerivRowIterator c_it = matrix.writeLine(cIndex++);
              c_it.addCol(i, typename CosseratNeedleSlidingConstraint<T>::Deriv(cy, vZero));
            }
            if (use[2]) {
              typename CosseratNeedleSlidingConstraint<T>::MatrixDerivRowIterator c_it = matrix.writeLine(cIndex++);
              c_it.addCol(i, typename CosseratNeedleSlidingConstraint<T>::Deriv(cz, vZero));
            }
        }
//        std::cout << "matrix : " << std::endl << matrix << std::endl;
        cMatrix.endEdit();
        self.m_nbLines = cIndex - self.m_constraintId;
    }


    template <class T>
    static void getConstraintViolation(CosseratNeedleSlidingConstraint<T>& self,
                                       const ConstraintParams* cParams,
                                       BaseVector *resV,
                                       const typename CosseratNeedleSlidingConstraint<T>::DataVecCoord &x,
                                       const typename CosseratNeedleSlidingConstraint<T>::DataVecDeriv &v)
    {
        if(self.d_componentState.getValue() != ComponentState::Valid)
            return ;

        SOFA_UNUSED(cParams);
        SOFA_UNUSED(x);
        SOFA_UNUSED(v);
        const typename CosseratNeedleSlidingConstraint<T>::VecCoord &positions = x.getValue();
        type::Vec<6,bool> use = self.d_useDirections.getValue();

        for (unsigned int i = 0; i < positions.size(); i++){
          if (use[1]) {
            typename CosseratNeedleSlidingConstraint<T>::Real dfree1 = positions[i][1];
            resV->set(self.m_constraintId + 2 * i, dfree1);
          }
          if (use[2]) {
            typename CosseratNeedleSlidingConstraint<T>::Real dfree2 = positions[i][2];
            resV->set(self.m_constraintId + 2 * i + 1, dfree2);
          }
        }
    }


    template<class T>
    static void getConstraintResolution(CosseratNeedleSlidingConstraint<T>& self,
                                        const ConstraintParams* cParams,
                                        std::vector<ConstraintResolution*>& resTab,
                                        unsigned int& offset)
    {
        SOFA_UNUSED(cParams);
        ReadAccessor<Data<typename CosseratNeedleSlidingConstraint<T>::VecCoord>> positions = self.mstate->readPositions();
        type::Vec<6,bool> use = self.d_useDirections.getValue();
        for (size_t i = 0; i < positions.size(); i++){
            if (use[1]) resTab[offset++] = new BilateralConstraintResolution();
            if (use[2]) resTab[offset++] = new BilateralConstraintResolution();
        }
    }

};


template <> SOFA_COSSERATPLUGIN_API
void CosseratNeedleSlidingConstraint<Rigid3Types>::buildConstraintMatrix(const ConstraintParams* cParams,
                                                                         DataMatrixDeriv &cMatrix,
                                                                         unsigned int &cIndex,
                                                                         const DataVecCoord &x)
{
    CosseratNeedleSlidingConstraintSpecialization<RigidImpl>::buildConstraintMatrix(*this, cParams, cMatrix, cIndex, x) ;
}


template <> SOFA_COSSERATPLUGIN_API
void CosseratNeedleSlidingConstraint<Rigid3Types>::getConstraintViolation(const ConstraintParams* cParams,
                                                                          BaseVector *resV,
                                                                          const DataVecCoord &x,
                                                                          const DataVecDeriv &v)
{
    CosseratNeedleSlidingConstraintSpecialization<RigidImpl>::getConstraintViolation(*this, cParams, resV, x, v);
}


template<> SOFA_COSSERATPLUGIN_API
void CosseratNeedleSlidingConstraint<Rigid3Types>::getConstraintResolution(const ConstraintParams* cParams,
                                                                           std::vector<ConstraintResolution*>& resTab,
                                                                           unsigned int& offset)
{
    CosseratNeedleSlidingConstraintSpecialization<RigidImpl>::getConstraintResolution(*this, cParams, resTab, offset);
}

int CosseratNeedleSlidingConstraintClass = RegisterObject("Simulate sliding contraints for needle insertion.")
    .add< CosseratNeedleSlidingConstraint<sofa::defaulttype::Vec3Types> >(true)
    .add< CosseratNeedleSlidingConstraint<sofa::defaulttype::Rigid3Types> >()
;

template class SOFA_COSSERATPLUGIN_API CosseratNeedleSlidingConstraint<Vec3Types>;
template class SOFA_COSSERATPLUGIN_API CosseratNeedleSlidingConstraint<Rigid3Types>;

} // namespace sofa


