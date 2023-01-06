#ifndef DrawTrianglesComponent_INL
#define DrawTrianglesComponent_INL

#include "DrawTrianglesComponent.h"
#include <stdexcept>
#include <exception>
#include <sofa/helper/AdvancedTimer.h>

#include <sofa/core/visual/VisualManager.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/FrameBufferObject.h>
#include <SofaOpenglVisual/OglShader.h>
#include <sofa/core/objectmodel/DataFileName.h>

#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseTopology/GridTopology.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/helper/decompose.h>
#include <sofa/helper/gl/template.h>
#include <assert.h>
#include <iostream>
#include <set>
#include <sofa/linearalgebra/CompressedRowSparseMatrix.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <cstdint>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>

namespace sofa::component::controller

DrawTrianglesComponent::DrawTrianglesComponent()
    : d_transparency(initData(&d_transparency, double(0.1), "transparency", "Scene scale")),
      d_vecTriangles(initData(&d_vecTriangles, "triangles", "tetra index elements")),
      d_vecTetra(initData(&d_vecTetra, "tetrahedra", "tetra index elements")),
      d_maxStress(initData(&d_maxStress, (double)5e3, "maxStress", "max Stress")),
      d_minStress(initData(&d_minStress, (double)0, "minStress", "min Stress")),
      d_maxVonMisesPerNode(initData(&d_maxVonMisesPerNode, "maxVonMisesPerNode", "max Von Mises Per Node")),
      d_draw(initData(&d_draw, true, "draw", "draw triangles"))
{
    this->f_listening.setValue(true);
}

void DrawTrianglesComponent::init()
{
    this->getContext()->get(m_state);
    this->getContext()->get(m_tetraForceField);
    if (m_state == NULL)
    {
        serr << "Error cannot find the mstate" << sendl;
        return;
    }
    if (m_tetraForceField == NULL)
    {
        serr << "Error cannot find the m_tetrahedron" << sendl;
        return;
    }
    m_VonMisesColorMap.setColorScheme("Blue to Red");
    m_VonMisesColorMap.reinit();
    m_tetraForceField->updateVonMisesStress = true;

    m_minVM = d_minStress.getValue();
    m_maxVM = d_maxStress.getValue();

    helper::ReadAccessor<Data<helper::vector<Real>>> vMN = m_tetraForceField->_vonMisesPerNode;
    helper::vector<double> &vonMises = *d_maxVonMisesPerNode.beginEdit();
    vonMises.resize(vMN.size());
    for (size_t i = 0; i < vMN.size(); i++)
    {
        vonMises[i] = vMN[i];
    }
    d_maxVonMisesPerNode.endEdit();
}

void DrawTrianglesComponent::handleEvent(sofa::core::objectmodel::Event *event)
{
    helper::ReadAccessor<Data<helper::vector<Real>>> vMN = m_tetraForceField->_vonMisesPerNode;
    if (dynamic_cast<sofa::simulation::AnimateBeginEvent *>(event))
    {
        helper::vector<double> &vonMises = *d_maxVonMisesPerNode.beginEdit();
        for (size_t index = 0; index < vMN.size(); index++)
        {
            m_minVM = (vMN[index] < m_minVM) ? vMN[index] : m_minVM;
            m_maxVM = (vMN[index] > m_maxVM) ? vMN[index] : m_maxVM;
            if (vonMises[index] < vMN[index])
            {
                vonMises[index] = vMN[index];
            }
        }
        d_maxVonMisesPerNode.endEdit();
    }
    else if (sofa::core::objectmodel::KeypressedEvent *ev = dynamic_cast<sofa::core::objectmodel::KeypressedEvent *>(event))
    {
        if ((ev->getKey() == 'l') || (ev->getKey() == 'L'))
        {
            helper::vector<double> &vonMises = *d_maxVonMisesPerNode.beginEdit();
            m_minVM = d_minStress.getValue();
            m_maxVM = d_maxStress.getValue();

            for (size_t i = 0; i < vMN.size(); i++)
            {
                vonMises[i] = 0.0;
            }
            d_maxVonMisesPerNode.endEdit();
        }
    }
}

void DrawTrianglesComponent::draw(const core::visual::VisualParams *vparams)
{
    if (d_draw.getValue())
        if (!this->m_state)
            return;
    m_tetraForceField->updateVonMisesStress = true;

    const VecCoord &x = this->m_state->read(core::ConstVecCoordId::position())->getValue();
    helper::ReadAccessor<Data<helper::vector<Real>>> vM = m_tetraForceField->_vonMisesPerElement;
    helper::ReadAccessor<Data<helper::vector<Real>>> vMN = m_tetraForceField->_vonMisesPerNode;
    // helper::vector<double> & vMN = *d_maxVonMisesPerNode.beginEdit();

    helper::ColorMap::evaluator<Real> evalColor = m_VonMisesColorMap.getEvaluator(d_minStress.getValue(), d_maxStress.getValue());

    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //    glColor4f(1,1,1,1);
    //    if (m_textureImage!=TEXTURE_UNASSIGNED) {
    //        glEnable(GL_TEXTURE_2D);
    //        glBindTexture(GL_TEXTURE_2D, m_textureImage);
    //        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    //    }
    if (vparams->displayFlags().getShowWireFrame())
    {
        glBegin(GL_LINES);

        int i = 0;
        std::vector<defaulttype::Vec4f> nodeColors(4);
        for (auto it = d_vecTriangles.getValue().begin(); it != d_vecTriangles.getValue().end(); ++it, ++i)
        {
            Index a = (*it)[0];
            Index b = (*it)[1];
            Index c = (*it)[2];
            nodeColors[0] = evalColor(vMN[a]);
            nodeColors[0][3] = d_transparency.getValue();
            nodeColors[1] = evalColor(vMN[b]);
            nodeColors[1][3] = d_transparency.getValue();
            nodeColors[2] = evalColor(vMN[c]);
            nodeColors[2][3] = d_transparency.getValue();

            float *ca = &nodeColors[0][0];
            float *cb = &nodeColors[1][0];
            float *cc = &nodeColors[2][0];

            Coord pa = x[a];
            Coord pb = x[b];
            Coord pc = x[c];

            glColor4fv(ca);
            helper::gl::glVertexT(pa);
            glColor4fv(cb);
            helper::gl::glVertexT(pb);
            glColor4fv(cc);
            helper::gl::glVertexT(pc);
        }
        glEnd();
    }
    else
    {

        glBegin(GL_TRIANGLES);

        //    for (size_t nd = 0; nd < x.size(); nd++) {
        //        defaulttype::Vec4f col = evalColor(vMN[nd]);

        //        glColor4f(col[0],col[1],col[2],1.0);
        //        glVertex3d(x[nd][0],x[nd][1],x[nd][2]);
        //    }

        // vparams->drawTool()->drawPoints(pts, 10, nodeColors);

        int i = 0;
        std::vector<defaulttype::Vec4f> nodeColors(4);

        i = 0;
        for (auto it = d_vecTetra.getValue().begin(); it != d_vecTetra.getValue().end(); ++it, ++i)
        {
            Index a = (*it)[0];
            Index b = (*it)[1];
            Index c = (*it)[2];
            Index d = (*it)[3];
            nodeColors[0] = evalColor(vMN[a]);
            nodeColors[0][3] = d_transparency.getValue();
            nodeColors[1] = evalColor(vMN[b]);
            nodeColors[1][3] = d_transparency.getValue();
            nodeColors[2] = evalColor(vMN[c]);
            nodeColors[2][3] = d_transparency.getValue();
            nodeColors[3] = evalColor(vMN[d]);
            nodeColors[3][3] = d_transparency.getValue();
            float *ca = &nodeColors[0][0];
            float *cb = &nodeColors[1][0];
            float *cc = &nodeColors[2][0];
            float *cd = &nodeColors[3][0];

            Coord center = (x[a] + x[b] + x[c] + x[d]) * 0.125;
            Coord pa = (x[a] + center) * (Real)0.666667;
            Coord pb = (x[b] + center) * (Real)0.666667;
            Coord pc = (x[c] + center) * (Real)0.666667;
            Coord pd = (x[d] + center) * (Real)0.666667;

            glColor4fv(ca);
            helper::gl::glVertexT(pa);
            glColor4fv(cb);
            helper::gl::glVertexT(pb);
            glColor4fv(cc);
            helper::gl::glVertexT(pc);

            glColor4fv(ca);
            helper::gl::glVertexT(pa);
            glColor4fv(cb);
            helper::gl::glVertexT(pb);
            glColor4fv(cd);
            helper::gl::glVertexT(pd);

            glColor4fv(ca);
            helper::gl::glVertexT(pa);
            glColor4fv(cc);
            helper::gl::glVertexT(pc);
            glColor4fv(cd);
            helper::gl::glVertexT(pd);

            glColor4fv(cc);
            helper::gl::glVertexT(pc);
            glColor4fv(cb);
            helper::gl::glVertexT(pb);
            glColor4fv(cd);
            helper::gl::glVertexT(pd);
        }

        for (auto it = d_vecTriangles.getValue().begin(); it != d_vecTriangles.getValue().end(); ++it, ++i)
        {
            Index a = (*it)[0];
            Index b = (*it)[1];
            Index c = (*it)[2];
            nodeColors[0] = evalColor(vMN[a]);
            nodeColors[0][3] = d_transparency.getValue();
            nodeColors[1] = evalColor(vMN[b]);
            nodeColors[1][3] = d_transparency.getValue();
            nodeColors[2] = evalColor(vMN[c]);
            nodeColors[2][3] = d_transparency.getValue();
            float *ca = &nodeColors[0][0];
            float *cb = &nodeColors[1][0];
            float *cc = &nodeColors[2][0];

            Coord pa = x[a];
            Coord pb = x[b];
            Coord pc = x[c];

            glColor4fv(ca);
            helper::gl::glVertexT(pa);
            glColor4fv(cb);
            helper::gl::glVertexT(pb);
            glColor4fv(cc);
            helper::gl::glVertexT(pc);
        }

        glEnd();
    }
}

} // end namespace controller

#endif
