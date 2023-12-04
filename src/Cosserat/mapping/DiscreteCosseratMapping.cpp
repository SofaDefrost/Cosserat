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
#define SOFA_COSSERAT_CPP_DiscreteCosseratMapping
#include "DiscreteCosseratMapping.inl"

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::component::mapping
{
using namespace sofa::defaulttype;


template <>
void DiscreteCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>:: applyJ(
    const core::MechanicalParams* /* mparams */, const type::vector< OutDataVecDeriv*>& dataVecOutVel,
    const type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
    const type::vector<const In2DataVecDeriv*>& dataVecIn2Vel) {

    if(d_debug.getValue())
        std::cout<< " ########## ApplyJ R Function ########"<< std::endl;

    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
        return;
    const In1VecDeriv& in1_vel = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& in2_vel = dataVecIn2Vel[0]->getValue();
    OutVecDeriv& out_vel = *dataVecOutVel[0]->beginEdit();
    const auto baseIndex = d_baseIndex.getValue();

    // Curv abscissa of nodes and frames
    helper::ReadAccessor<Data<List>> curv_abs_section =  d_curv_abs_section;
    helper::ReadAccessor<Data<List>> curv_abs_frames = d_curv_abs_frames;

    const In1VecCoord& inDeform = m_fromModel1->read(core::ConstVecCoordId::position())->getValue(); //Strains
    // Compute the tangent Exponential SE3 vectors
    this->update_TangExpSE3(inDeform);

    //Get base velocity as input this is also called eta
    m_nodesVelocityVectors.clear();

    //Get base velocity and convert to Vec6, for the facility of computation
    type::Vec6 baseVelocity; //
    for (auto u=0; u<6; u++)
        baseVelocity[u] = in2_vel[baseIndex][u];

    //Apply the local transform i.e. from SOFA's frame to Cosserat's frame
    const In2VecCoord& xfrom2Data = m_fromModel2->read(core::ConstVecCoordId::position())->getValue();
    Transform TInverse = Transform(xfrom2Data[baseIndex].getCenter(), xfrom2Data[baseIndex].getOrientation()).inversed();
    Mat6x6 P = this->build_projector(TInverse);
    type::Vec6 baseLocalVelocity = P * baseVelocity; //This is the base velocity in Locale frame
    m_nodesVelocityVectors.push_back(baseLocalVelocity);
    if(d_debug.getValue())
        std::cout << "Base local Velocity :"<< baseLocalVelocity <<std::endl;

    //Compute velocity at nodes
    for (unsigned int i = 1 ; i < curv_abs_section.size(); i++) {
        Transform Trans = m_nodesExponentialSE3Vectors[i].inversed();
        Tangent Adjoint; Adjoint.clear();
        this->computeAdjoint(Trans, Adjoint);

        type::Vec6 node_Xi_dot;
        for (unsigned int u =0; u<6; u++)
            node_Xi_dot(i) = in1_vel[i-1][u];

        Vec6 eta_node_i = Adjoint * (m_nodesVelocityVectors[i-1] + m_nodesTangExpVectors[i] *node_Xi_dot );
        m_nodesVelocityVectors.push_back(eta_node_i);
        if(d_debug.getValue())
            std::cout<< "Node velocity : "<< i << " = " << eta_node_i<< std::endl;
    }
    const OutVecCoord& out = m_toModel->read(core::ConstVecCoordId::position())->getValue();
    auto sz = curv_abs_frames.size();
    out_vel.resize(sz);

    for (unsigned int i = 0 ; i < sz; i++) {
        Transform Trans = m_framesExponentialSE3Vectors[i].inversed();
        Tangent Adjoint; Adjoint.clear();
        this->computeAdjoint(Trans, Adjoint);
        type::Vec6 frame_Xi_dot = in1_vel[m_indicesVectors[i]-1];
//        for (unsigned int u =0; u<6; u++)
//            frame_Xi_dot(i) = in1_vel[m_indicesVectors[i]-1][u];

        Vec6 eta_frame_i = Adjoint * (m_nodesVelocityVectors[m_indicesVectors[i]-1] + m_framesTangExpVectors[i] * frame_Xi_dot ); // eta

        auto T = Transform(out[i].getCenter(), out[i].getOrientation());
        Tangent Proj = this->build_projector(T);

        out_vel[i] = Proj * eta_frame_i;
        if(d_debug.getValue())
            std::cout<< "Frame velocity : "<< i << " = " << eta_frame_i<< std::endl;
    }
    dataVecOutVel[0]->endEdit();
    m_index_input = 0;
}



template <>
void DiscreteCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>:: applyJT(
    const core::MechanicalParams* /*mparams*/, const type::vector< In1DataVecDeriv*>& dataVecOut1Force,
    const type::vector< In2DataVecDeriv*>& dataVecOut2Force,
    const type::vector<const OutDataVecDeriv*>& dataVecInForce)  {

    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

    if(d_debug.getValue())
        std::cout<< " ########## ApplyJT force R Function ########"<< std::endl;
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

    Tangent matB_trans; matB_trans.clear();
    for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;


    for (auto s = sz ; s-- ; ) {
        Tangent coAdjoint;

        this->compute_coAdjoint(m_framesExponentialSE3Vectors[s], coAdjoint);  // m_framesExponentialSE3Vectors[s] computed in apply
        Vec6 node_F_Vec = coAdjoint * local_F_Vec[s];
        Mat6x6 temp = m_framesTangExpVectors[s];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
        temp.transpose();
        type::Vec6 f = matB_trans * temp * node_F_Vec;

        if(index != m_indicesVectors[s]){
            index--;
            //bring F_tot to the reference of the new beam
            this->compute_coAdjoint(m_nodesExponentialSE3Vectors[index],coAdjoint);  //m_nodesExponentialSE3Vectors computed in apply
            F_tot = coAdjoint * F_tot;
            Mat6x6 temp = m_nodesTangExpVectors[index];
            temp.transpose();
            //apply F_tot to the new beam
            type::Vec6 temp_f = matB_trans * temp * F_tot;
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

template <>
void DiscreteCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>::applyJT(
    const core::ConstraintParams*/*cparams*/ , const type::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
    const type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
    const type::vector<const OutDataMatrixDeriv*>& dataMatInConst)
{
    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
        return;

    if(d_debug.getValue())
        std::cout<< " ########## ApplyJT constraint R Function ########"<< std::endl;

    //We need only one input In model and input Root model (if present)
    In1MatrixDeriv& out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the strain space (reduced coordinate)
    In2MatrixDeriv& out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the reference frame (base frame)
    const OutMatrixDeriv& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames

    const OutVecCoord& frame = m_toModel->read(core::ConstVecCoordId::position())->getValue();
    const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
    const In1VecCoord x1from = x1fromData->getValue();

    Tangent matB_trans; matB_trans.clear();
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
        auto colIt = rowIt.begin();
        auto colItEnd = rowIt.end();

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

            type::Vec6 valueConst;
            for(unsigned j = 0; j < 6; j++)
                valueConst[j] = colIt.val()[j];

            int indexBeam =  m_indicesVectors[childIndex];

            Transform _T = Transform(frame[childIndex].getCenter(),frame[childIndex].getOrientation());
            Tangent P_trans =(this->build_projector(_T));
            P_trans.transpose();

            Mat6x6 coAdjoint;
            this->compute_coAdjoint(m_framesExponentialSE3Vectors[childIndex], coAdjoint);  // m_framesExponentialSE3Vectors[s] computed in apply
            Tangent temp = m_framesTangExpVectors[childIndex];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
            temp.transpose();

            type::Vec6 local_F =  coAdjoint * P_trans * valueConst; // constraint direction in local frame of the beam.

            type::Vec6 f = matB_trans * temp * local_F; // constraint direction in the strain space.

            o1.addCol(indexBeam-1, f);
            std::tuple<int,Vec6> node_force = std::make_tuple(indexBeam, local_F);

            NodesInvolved.push_back(node_force);
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
                type::Vec6 temp_f = matB_trans * temp * CumulativeF;

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


// Register in the Factory
int DiscreteCosseratMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a rigid parent")
                                       .add< DiscreteCosseratMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types > >(true)
                                       .add< DiscreteCosseratMapping< sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types > >()
    ;

template class SOFA_COSSERAT_API DiscreteCosseratMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
template class SOFA_COSSERAT_API DiscreteCosseratMapping< sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;

} // namespace sofa::component::mapping
