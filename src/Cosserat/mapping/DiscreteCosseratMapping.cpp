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
#include <Cosserat/mapping/DiscreteCosseratMapping.inl>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace Cosserat::mapping
{
using namespace sofa::defaulttype;


template <>
void DiscreteCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>:: applyJ(
    const sofa::core::MechanicalParams* /* mparams */, const vector< sofa::DataVecDeriv_t<Out>*>& dataVecOutVel,
    const vector<const sofa::DataVecDeriv_t<In1>*>& dataVecIn1Vel,
    const vector<const sofa::DataVecDeriv_t<In2>*>& dataVecIn2Vel) {

    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
        return;

    msg_info() << " ########## ApplyJ R Function ########";

    const sofa::VecDeriv_t<In1>& in1_vel = sofa::helper::getReadAccessor(*dataVecIn1Vel[0]);
    const sofa::VecDeriv_t<In2>& in2_vel = sofa::helper::getReadAccessor(*dataVecIn2Vel[0]);

    auto out_vel = sofa::helper::getWriteOnlyAccessor(*dataVecOutVel[0]);

    const auto baseIndex = d_baseIndex.getValue();

    // Curv abscissa of nodes and frames
    const auto curv_abs_section = sofa::helper::getReadAccessor(d_curv_abs_section);
    const auto curv_abs_frames = sofa::helper::getReadAccessor(d_curv_abs_frames);

    const auto inDeform = sofa::helper::getReadAccessor(*m_strain_state->read(sofa::core::vec_id::read_access::position));

    // Compute the tangent Exponential SE3 vectors
    this->updateTangExpSE3(inDeform);

    //Get base velocity as input this is also called eta
    m_nodesVelocityVectors.clear();

    //Get base velocity and convert to Vec6, for the facility of computation
    Vec6 baseVelocity; //
    for (auto u=0; u<6; u++)
        baseVelocity[u] = in2_vel[baseIndex][u];

    //Apply the local transform i.e. from SOFA's frame to Cosserat's frame
    const sofa::VecCoord_t<In2>& xfrom2Data = sofa::helper::getReadAccessor(*m_rigid_base->read(sofa::core::vec_id::read_access::position));
    Frame TInverse = Frame(xfrom2Data[baseIndex].getCenter(), xfrom2Data[baseIndex].getOrientation()).inversed();
    Mat6x6 P = this->buildProjector(TInverse);
    Vec6 baseLocalVelocity = P * baseVelocity; //This is the base velocity in Locale frame
    m_nodesVelocityVectors.push_back(baseLocalVelocity);

    msg_info() << "Base local Velocity :"<< baseLocalVelocity;

    //Compute velocity at nodes
    for (unsigned int i = 1 ; i < curv_abs_section.size(); i++)
    {
        Frame Trans = m_nodesExponentialSE3Vectors[i].inversed();
        TangentTransform Adjoint;
        this->computeAdjoint(Trans, Adjoint);

        Vec6 node_Xi_dot;
        for (unsigned int u =0; u<6; u++)
            node_Xi_dot(i) = in1_vel[i-1][u];

        Vec6 eta_node_i = Adjoint * (m_nodesVelocityVectors[i-1] + m_nodesTangExpVectors[i] *node_Xi_dot );
        m_nodesVelocityVectors.push_back(eta_node_i);
        msg_info() << "Node velocity : "<< i << " = " << eta_node_i;
    }
    const sofa::VecCoord_t<Out>& out = sofa::helper::getReadAccessor(*m_global_frames->read(sofa::core::vec_id::read_access::position));

    auto sz = curv_abs_frames.size();
    out_vel.resize(sz);
    for (unsigned int i = 0 ; i < sz; i++) {
        Frame Trans = m_framesExponentialSE3Vectors[i].inversed();
        TangentTransform Adjoint; Adjoint.clear();
        this->computeAdjoint(Trans, Adjoint);
        Vec6 frame_Xi_dot = in1_vel[m_indicesVectors[i]-1];
        Vec6 eta_frame_i = Adjoint * (m_nodesVelocityVectors[m_indicesVectors[i]-1] + m_framesTangExpVectors[i] * frame_Xi_dot ); // eta

        auto T = Frame(out[i].getCenter(), out[i].getOrientation());
        TangentTransform Proj = this->buildProjector(T);

        out_vel[i] = Proj * eta_frame_i;
        msg_info() << "Frame velocity : "<< i << " = " << eta_frame_i;
    }
    m_indexInput = 0;
}

template <>
void DiscreteCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>:: applyJT(
    const sofa::core::MechanicalParams* /*mparams*/, const vector< sofa::DataVecDeriv_t<In1>*>& dataVecOut1Force,
    const vector< sofa::DataVecDeriv_t<In2>*>& dataVecOut2Force,
    const vector<const sofa::DataVecDeriv_t<Out>*>& dataVecInForce)  {

    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

    msg_info() << " ########## ApplyJT force R Function ########";
    const sofa::VecDeriv_t<Out>& in = dataVecInForce[0]->getValue();

    sofa::VecDeriv_t<In1> out1 = sofa::helper::getWriteAccessor(*dataVecOut1Force[0]);
    sofa::VecDeriv_t<In2> out2 = sofa::helper::getWriteAccessor(*dataVecOut2Force[0]);
    const auto baseIndex = d_baseIndex.getValue();

    const sofa::VecCoord_t<Out>& frame =
        m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();
    const sofa::DataVecCoord_t<In1>* x1fromData =
        m_strain_state->read(sofa::core::vec_id::read_access::position);
    const sofa::VecCoord_t<In1> x1from = x1fromData->getValue();
    vector<Vec6> local_F_Vec;  local_F_Vec.clear();

    out1.resize(x1from.size());

    //convert the input from Deriv type to vec6 type, for the purpose of the matrix vector multiplication
    for (unsigned int var = 0; var < in.size(); ++var) {
        Vec6 vec;
        for(unsigned j = 0; j < 6; j++) vec[j] = in[var][j];
        //Convert input from global frame(SOFA) to local frame
        Frame _T = Frame(frame[var].getCenter(),frame[var].getOrientation());
        Mat6x6 P_trans =(this->buildProjector(_T)); P_trans.transpose();
        Vec6 local_F = P_trans * vec;
        local_F_Vec.push_back(local_F);
    }

    //Compute output forces
    auto sz = m_indicesVectors.size();
    auto index =  m_indicesVectors[sz-1];
    m_totalBeamForceVectors.clear();
    m_totalBeamForceVectors.resize(sz);

    Vec6 F_tot; F_tot.clear();
    m_totalBeamForceVectors.push_back(F_tot);

    TangentTransform matB_trans; matB_trans.clear();
    for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;

    for (auto s = sz ; s-- ; ) {
        TangentTransform coAdjoint;

        this->computeCoAdjoint(m_framesExponentialSE3Vectors[s], coAdjoint);  // m_framesExponentialSE3Vectors[s] computed in apply
        Vec6 node_F_Vec = coAdjoint * local_F_Vec[s];
        Mat6x6 temp = m_framesTangExpVectors[s];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
        temp.transpose();
        Vec6 f = matB_trans * temp * node_F_Vec;

        if(index != m_indicesVectors[s]){
            index--;
            //bring F_tot to the reference of the new beam
            this->computeCoAdjoint(m_nodesExponentialSE3Vectors[index],coAdjoint);  //m_nodesExponentialSE3Vectors computed in apply
            F_tot = coAdjoint * F_tot;
            Mat6x6 temp_mat = m_nodesTangExpVectors[index];
            temp_mat.transpose();
            //apply F_tot to the new beam
            const Vec6 temp_f = matB_trans * temp_mat * F_tot;
            out1[index-1] += temp_f;
        }

        msg_info() << "f at s ="<< s <<" and index"<< index <<  " is : "<< f;

        //compute F_tot
        F_tot += node_F_Vec;
        out1[m_indicesVectors[s]-1] += f;
    }

    const auto frame0 = Frame(frame[0].getCenter(),frame[0].getOrientation());
    const Mat6x6 M = this->buildProjector(frame0);
    out2[baseIndex] += M * F_tot;

    msg_info()
            << "Node forces "<< out1 << msgendl
            << "base Force: "<< out2[baseIndex];
}

template <>
void DiscreteCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>::applyJT(
    const sofa::core::ConstraintParams*/*cparams*/ , const vector< sofa::DataMatrixDeriv_t<In1>*>&  dataMatOut1Const,
    const vector< sofa::DataMatrixDeriv_t<In2>*>&  dataMatOut2Const ,
    const vector<const sofa::DataMatrixDeriv_t<Out>*>& dataMatInConst)
{
    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
        return;

    msg_info() << " ########## ApplyJT constraint R Function ########";

    //We need only one input In model and input Root model (if present)
    sofa::MatrixDeriv_t<In1>& out1 = sofa::helper::getWriteAccessor(*dataMatOut1Const[0]); // constraints on the strain space (reduced coordinate)
    sofa::MatrixDeriv_t<In2>& out2 = sofa::helper::getWriteAccessor(*dataMatOut2Const[0]); // constraints on the reference frame (base frame)
    const sofa::MatrixDeriv_t<Out>& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames

    const sofa::VecCoord_t<Out>& frame =
        m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();
    const sofa::DataVecCoord_t<In1>* x1fromData =
        m_strain_state->read(sofa::core::vec_id::read_access::position);
    const sofa::VecCoord_t<In1> x1from = x1fromData->getValue();

    TangentTransform matB_trans; matB_trans.clear();
    for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;

    vector< std::tuple<int,Vec6> > NodesInvolved;
    vector< std::tuple<int,Vec6> > NodesInvolvedCompressed;

    sofa::MatrixDeriv_t<Out>::RowConstIterator rowItEnd = in.end();

    bool doPrintLog = f_printLog.getValue();
    for (sofa::MatrixDeriv_t<Out>::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        msg_info_when(doPrintLog)
             <<"************* Apply JT (MatrixDeriv) iteration on line "
             << rowIt.index()
             <<"*************  ";

        auto colIt = rowIt.begin();
        auto colItEnd = rowIt.end();

        // Creates a constraints if the input constraint is not empty.
        if (colIt == colItEnd)
        {
            msg_info_when(doPrintLog)
                <<"no column for this constraint";
            continue;
        }
        sofa::MatrixDeriv_t<In1>::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
        sofa::MatrixDeriv_t<In2>::RowIterator o2 = out2.writeLine(rowIt.index());

        NodesInvolved.clear();
        while (colIt != colItEnd)
        {
            int childIndex = colIt.index();

            Vec6 valueConst;
            for(unsigned j = 0; j < 6; j++)
                valueConst[j] = colIt.val()[j];

            int indexBeam =  m_indicesVectors[childIndex];

            const auto _T = Frame(frame[childIndex].getCenter(),frame[childIndex].getOrientation());
            TangentTransform P_trans =(this->buildProjector(_T));
            P_trans.transpose();

            Mat6x6 coAdjoint;
            this->computeCoAdjoint(m_framesExponentialSE3Vectors[childIndex], coAdjoint);  // m_framesExponentialSE3Vectors[s] computed in apply
            Mat6x6 temp = m_framesTangExpVectors[childIndex];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
            temp.transpose();

            Vec6 local_F =  coAdjoint * P_trans * valueConst; // constraint direction in local frame of the beam.

            Vec6 f = matB_trans * temp * local_F; // constraint direction in the strain space.

            o1.addCol(indexBeam-1, f);
            std::tuple<int,Vec6> node_force = std::make_tuple(indexBeam, local_F);

            NodesInvolved.push_back(node_force);
            colIt++;
        }

        if(doPrintLog)
        {
            std::stringstream tmp;
            for (size_t i = 0; i < NodesInvolved.size(); i++)
                tmp << "index :" <<get<0>(NodesInvolved[i]) << " force :" << get<1>(NodesInvolved[i]) << msgendl;
            msg_info() <<"==> NodesInvolved : " << tmp.str();
        }

        // sort the Nodes Involved by decreasing order
        std::sort(begin(NodesInvolved), end(NodesInvolved),
                  [](std::tuple<int, Vec6> const &t1, std::tuple<int, Vec6> const &t2) {
                      return std::get<0>(t1) > std::get<0>(t2);
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

        if(doPrintLog)
        {
            std::stringstream tmp;
            tmp<<" NodesInvolved after sort and compress"<<msgendl;
            for (size_t i = 0; i < NodesInvolvedCompressed.size(); i++)
                tmp << "index :" <<get<0>(NodesInvolvedCompressed[i]) << " force :"
                          << get<1>(NodesInvolvedCompressed[i]) << msgendl;
            msg_info() << tmp.str();
        }

        auto baseIndex = d_baseIndex.getValue();
        for (unsigned n=0; n<NodesInvolvedCompressed.size(); n++)
        {
            std::tuple<int,Vec6> test = NodesInvolvedCompressed[n];
            int numNode= std::get<0>(test);
            int i = numNode;
            Vec6 CumulativeF = std::get<1>(test);

            while(i>0)
            {
                //cumulate on beam frame
                Mat6x6 co_adjoint;
                this->computeCoAdjoint(m_nodesExponentialSE3Vectors[i-1],co_adjoint);
                CumulativeF = co_adjoint * CumulativeF;
                // transfer to strain space (local coordinates)
                Mat6x6 temp = m_nodesTangExpVectors[i-1];
                temp.transpose();
                Vec6 temp_f = matB_trans * temp * CumulativeF;

                if(i>1) o1.addCol(i-2, temp_f);
                i--;
            }

            const auto frame0 = Frame(frame[0].getCenter(),frame[0].getOrientation());
            const Mat6x6 M = this->buildProjector(frame0);

            Vec6 base_force = M * CumulativeF;
            o2.addCol(baseIndex, base_force);
        }
    }
}


template class SOFA_COSSERAT_API DiscreteCosseratMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
template class SOFA_COSSERAT_API DiscreteCosseratMapping< sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;

} // namespace sofa::component::mapping

namespace Cosserat
{
// Register in the Factory
void registerDiscreteCosseratMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData(
        "This component facilitates the creation of Cosserat Cables in SOFA simulations. It takes two mechanical"
        "objects as inputs: the rigid base of the beam (with 6 degrees of freedom) and the local coordinates of the beam. Using "
        "these inputs, the component computes and outputs the mechanical positions of the beam in global coordinates. "
        "Like any mapping, it updates the positions and velocities of the outputs based on the inputs. "
        "Additionally, forces applied to the outputs are propagated back to the inputs, ensuring bidirectional coupling.")
        .add< mapping::DiscreteCosseratMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types > >(true)
        .add< mapping::DiscreteCosseratMapping< sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types > >());
}
}