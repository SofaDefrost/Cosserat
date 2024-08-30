# -*- coding: utf-8 -*-
"""
Cosserat class in SofaPython3.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

import Sofa
from useful.utils import addEdgeCollision, addPointsCollision
from useful.header import addHeader, addVisual, addSolverNode
from useful.params import Parameters, BeamGeometryParameters
from useful.geometry import CosseratGeometry, generate_edge_list
from numpy import array


class CosseratBase(Sofa.Prefab):
    """
    CosseratBase model prefab class. It is a prefab class that allow to create a cosserat beam/rod in Sofa.
           Structure:
           Node : {
                name : 'CosseratBase'
                Node0 MechanicalObject :     // Rigid position of the base of the beam
                Node1 MechanicalObject :    // Vec3d, cosserat local parameters composed of the twist and the bending along y and z
                Node1 ForceField          // Base on Hook's law, it computed the force applied on the beam
                (Node0-Node1)-child MechanicalObject     //  The child of the two precedent nodes, Rigid positions
                Allow to compute the cosserat frame in the world frame (Sofa frame)
                    Cosserat Mapping  //  it allow the transfer from the locial to the word frame
            }
            params

    """

    prefabParameters = [
        {"name": "name", "type": "string", "help": "Node name", "default": "Cosserat"},
        {
            "name": "translation",
            "type": "Vec3d",
            "help": "Cosserat base Rotation",
            "default": array([0.0, 0.0, 0.0]),
        },
        {
            "name": "rotation",
            "type": "Vec3d",
            "help": "Cosserat base Rotation",
            "default": array([0.0, 0.0, 0.0]),
        },
        {  # @TODO: to be removed
            "name": "attachingToLink",
            "type": "string",
            "help": "a rest shape force field will constraint the object "
            "to follow arm position",
            "default": "1",
        },
        {
            "name": "showObject",
            "type": "string",
            "help": " Draw object arrow ",
            "default": "0",
        },
    ]

    def __init__(self, *args, **kwargs):
        Sofa.Prefab.__init__(self, *args, **kwargs)
        self.params = kwargs.get(
            "params", Parameters()
        )  # Use the Parameters class with default values
        beamPhysicsParams = self.params.beamPhysicsParams
        beamGeometryParams = self.params.beamGeoParams

        self.beamMass = beamPhysicsParams.beamMass  # self.cosseratGeometry['beamMass']
        # Todo: add option in case None
        self.parent = kwargs.get("parent")
        self.useInertiaParams = beamPhysicsParams.useInertia  # False
        self.radius = beamPhysicsParams.beamRadius  # kwargs.get('radius')

        # @TODO: To be removed
        if self.parent.hasObject("EulerImplicitSolver") is False:
            print("The code does not have parent EulerImplicit")
            self.solverNode = addSolverNode(self.parent)
        else:
            self.solverNode = self.parent

        if "inertialParams" in kwargs:
            self.useInertiaParams = True
            self.inertialParams = kwargs["inertialParams"]

        self.rigidBaseNode = self._addRigidBaseNode()

        cosserat_geometry = CosseratGeometry(beamGeometryParams)
        self.frames3D = cosserat_geometry.cable_positionF

        self.cosseratCoordinateNode = self._addCosseratCoordinate(
            cosserat_geometry.bendingState, cosserat_geometry.sectionsLengthList
        )
        self.cosseratFrame = self._addCosseratFrame(
            cosserat_geometry.framesF,
            cosserat_geometry.curv_abs_inputS,
            cosserat_geometry.curv_abs_outputF,
        )

    def init(self):
        pass

    def addCollisionModel(self):
        tab_edges = generate_edge_list(self.frames3D)
        return addEdgeCollision(self.cosseratFrame, self.frames3D, tab_edges)

    def _addPointCollisionModel(self, nodeName="CollisionPoints"):
        tab_edges = generate_edge_list(self.frames3D)
        return addPointsCollision(
            self.cosseratFrame, self.frames3D, tab_edges, nodeName
        )

    def _addSlidingPoints(self):
        slidingPoint = self.cosseratFrame.addChild("slidingPoint")
        slidingPoint.addObject("MechanicalObject", name="slidingPointMO", position=self.frames3D, 
                               showObject="0", showIndices="0")
        slidingPoint.addObject("IdentityMapping")
        return slidingPoint

    def _addSlidingPointsWithContainer(self):
        slidingPoint = self.cosseratFrame.addChild("slidingPoint")
        slidingPoint.addObject("PointSetTopologyContainer")
        slidingPoint.addObject("PointSetTopologyModifier")
        slidingPoint.addObject(
            "MechanicalObject",
            name="slidingPointMO",
            position=self.frames3D,
            showObject="1",
            showIndices="0",
        )
        slidingPoint.addObject("IdentityMapping")
        return slidingPoint

    def _addRigidBaseNode(self):
        rigidBaseNode = self.addChild("rigidBase")
        # To be improved with classes in top
        positions = [[self.params.beamGeoParams.init_pos] + [0.0, 0.0, 0.0, 1.0]]

        rigidBaseNodeMo = rigidBaseNode.addObject(
            "MechanicalObject",
            template="Rigid3d",
            name="RigidBaseMO",
            showObjectScale=0.2,
            position=positions,
            translation=self.translation,
            rotation=self.rotation
        )
        rigidBaseNodeMo.showObject.setParent(self.showObject)

        # @TODO: remove this hard coded.
        # one can choose to set this to false and directly attach the beam base
        # to a control object in order to be able to drive it.
        if int(self.attachingToLink.value):
            print("Adding the rest shape to the base")
            rigidBaseNode.addObject(
                "RestShapeSpringsForceField",
                name="spring",
                stiffness=1e8,
                angularStiffness=1.0e8,
                external_points=0,
                mstate="@RigidBaseMO",
                points=0,
                template="Rigid3d",
            )
        return rigidBaseNode

    def _addCosseratCoordinate(self, bendingStates, listOfSectionsLength):
        cosseratCoordinateNode = self.addChild("cosseratCoordinate")
        cosseratCoordinateNode.addObject(
            "MechanicalObject",
            template="Vec3d",
            name="cosseratCoordinateMO",
            position=bendingStates,
            showIndices=0,
        )

        if self.useInertiaParams is False:
            cosseratCoordinateNode.addObject(
                "BeamHookeLawForceField",
                crossSectionShape=self.params.beamPhysicsParams.beamShape,
                length=listOfSectionsLength,
                radius=self.radius,
                youngModulus=self.params.beamPhysicsParams.youngModulus,
                poissonRatio=self.params.beamPhysicsParams.poissonRatio,
                rayleighStiffness=self.params.simuParams.rayleighStiffness,
                lengthY=self.params.beamPhysicsParams.length_Y,
                lengthZ=self.params.beamPhysicsParams.length_Z,
            )
        else:
            self._extracted_from_addCosseratCoordinate_15(
                cosseratCoordinateNode, listOfSectionsLength
            )
        return cosseratCoordinateNode

    # TODO Rename this here and in `addCosseratCoordinate`
    def _extracted_from_addCosseratCoordinate_15(
        self, cosseratCoordinateNode, listOfSectionsLength
    ):
        GA = self.params.beamPhysicsParams.GA
        GI = self.params.beamPhysicsParams.GI
        EA = self.params.beamPhysicsParams.EA
        EI = self.params.beamPhysicsParams.EI
        cosseratCoordinateNode.addObject(
            "BeamHookeLawForceField",
            crossSectionShape=self.params.beamPhysicsParams.beamShape,
            length=listOfSectionsLength,
            radius=self.params.beamPhysicsParams.beamRadius,
            useInertiaParams=True,
            GI=GI,
            GA=GA,
            EI=EI,
            EA=EA,
            rayleighStiffness=self.rayleighStiffness.value,
            lengthY=self.params.beamPhysicsParams.length_Y,
            lengthZ=self.params.beamPhysicsParams.length_Z,
        )

    def _addCosseratFrame(self, framesF, curv_abs_inputS, curv_abs_outputF):
        cosseratInSofaFrameNode = self.rigidBaseNode.addChild("cosseratInSofaFrameNode")
        self.cosseratCoordinateNode.addChild(cosseratInSofaFrameNode)
        framesMO = cosseratInSofaFrameNode.addObject(
            "MechanicalObject",
            template="Rigid3d",
            name="FramesMO",
            position=framesF,
            showIndices=self.params.beamGeoParams.show_frames_indices,
            showObject=self.params.beamGeoParams.show_frames_object,
            showObjectScale=1.8,  # Todo: remove this hard code
        )

        cosseratInSofaFrameNode.addObject(
            "UniformMass", totalMass=self.beamMass, showAxisSizeFactor="0"
        )

        cosseratInSofaFrameNode.addObject(
            "DiscreteCosseratMapping",
            curv_abs_input=curv_abs_inputS,
            curv_abs_output=curv_abs_outputF,
            name="cosseratMapping",
            input1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
            input2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
            output=framesMO.getLinkPath(),
            debug=0,
            radius=self.radius,
        )
        return cosseratInSofaFrameNode


Params = Parameters(beamGeoParams=BeamGeometryParameters(init_pos=[0, 0, 0]))


def createScene(rootNode):
    addHeader(rootNode)
    addVisual(rootNode)

    rootNode.findData("dt").value = 0.01
    rootNode.findData("gravity").value = [0.0, -9.81, 0.0]
    rootNode.addObject("BackgroundSetting", color="0 0.168627 0.211765")
    rootNode.addObject("FreeMotionAnimationLoop")
    rootNode.addObject("GenericConstraintSolver", tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject("Camera", position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild("solverNode")
    solverNode.addObject(
        "EulerImplicitSolver", rayleighStiffness="0.2", rayleighMass="0.1"
    )
    solverNode.addObject(
        "SparseLDLSolver", name="solver", template="CompressedRowSparseMatrixd"
    )
    solverNode.addObject("GenericConstraintCorrection")

    # Create a
    cosserat = solverNode.addChild(CosseratBase(parent=solverNode, params=Params))

    # use this to add the collision if the beam will interact with another object
    cosserat.addCollisionModel()

    # Attach a force at the beam tip,
    # we can attach this force to a non-mechanical node to control the beam in order to be able to drive it.
    cosserat.cosseratFrame

    return rootNode
