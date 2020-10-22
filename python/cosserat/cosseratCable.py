from math import sqrt, pi
from splib.objectmodel import SofaPrefab, SofaObject
from splib.numerics import Vec3, Quat
from splib.animation import animate, AnimationManager
from stlib.physics.deformable import ElasticMaterialObject
from stlib.physics.constraints import FixedBox
from stlib.scene import Node
#from cosseratUtilities import compute_BeamLenght, createCurvAbsOutput, createFramesList, extractFEMConstraintPoints
import Sofa
import os
from grippercontroller import GripperController
path = os.path.dirname(os.path.abspath(__file__))+'/../../scenes/inverseModelScenes/mesh/'

FEMpos = [" 0. 0. 0. 15. 0. 0. 30. 0. 0. 45. 0. 0. 60. 0. 0. 66. 0. 0. 81. 0.0 0.0"]

def cable(
    attachedTo=None,
    name="rigidBase",
    cableLenght= 81, 
    numberBeams = 6,
    stiffness="50000",
    angularStiffness=50000,
    rotation=[0.0,0.0,0.0,1.0],
    translation=[0.0,0.0,0.0],
    listOfBeamslenght = None,
    crossSectionShape='circular',
    radius=0.5,
    youngModulus=5e6    
    ):
  
    #  This create a new node in the scene. This node is appended to the finger's node.
    cable = attachedTo.createChild(name)
    
    # ###############
    # RigidBase
    ###############
    base = [translation[0],translation[1],translation[2],rotation[0],rotation[1],rotation[2],rotation[3]]
    RigidBaseMO = cable.createObject('MechanicalObject', template='Rigid3d',name="RigidBaseMO", position=base, showObject='1', showObjectScale='0.1')
    cable.createObject('RestShapeSpringsForceField', name='spring', stiffness=stiffness,angularStiffness=angularStiffness, external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
    
    
    ###############
    # Rate of angular Deformation  (2 sections)
    ###############
    #ratePosition = []
    #lenghts = []
    #if listOfBeamslenght != None:
        #lenghts = listOfBeamslenght
    #else:
        #beamLenght = cableLenght/numberBeams 
        #for i in range(0,numberBeams):
            #lenghts.append([beamLenght])
            
    #for i in range(0,numberBeams):
        #ratePosition.append([0.,0.,0.])
        
    #rateAngularDeformNode = attachedTo.createChild('rateAngularDeform')
    #rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=ratePosition)
    #BeamHookeLawForce = rateAngularDeformNode.createObject('BeamHookeLawForceField', crossSectionShape=crossSectionShape, length=lenghts, radius=radius, youngModulus=youngModulus)
        
    #return cable


    
    #mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, input2=inputCableMO, output=outputPointMO, direction=direction)

def addConstraintPoints(attachedTo, cstPoints,mappedPointsNode,translation=[0.,0.,0.],rotation=[0.,0.,0.]):
        
        trunkMappedPoints = attachedTo.createChild('constraintPoints')        
        inputFEMCable = trunkMappedPoints.createObject('MechanicalObject', name="pointsInFEM", position=cstPoints, 
                                                       showObject="1", showIndices="1", translation=translation, rotation=rotation)
        
        trunkMappedPoints.addChild(mappedPointsNode)
        trunkMappedPoints.createObject('BarycentricMapping')
        return inputFEMCable.getLinkPath()


def Finger(parentNode=None, name="Finger",
           rotation=[0.0, 0.0, 0.0], translation=[0.0, 0.0, 0.0],
           fixingBox=[-18., -15., -8., 2., -3., 8]):
    
    finger = ElasticMaterialObject(parentNode,name=name,
                                   volumeMeshFileName=path+"transFinger.vtk",
                                   poissonRatio=0.45,
                                   youngModulus=600,
                                   totalMass=0.5,
                                   surfaceColor=[0.0, 0.7, 0.7],
                                   surfaceMeshFileName=path+"transFinger.stl",
                                   rotation=rotation,
                                   translation=translation)

    FixedBox(finger.node, atPositions=fixingBox, doVisualization=True)
    
    #CollisionMesh(eobject.node, name="CollisionMesh",
                 #surfaceMeshFileName=path+"finger.stl",
                 #rotation=rotation, translation=translation,
                 #collisionGroup=[1, 2])

    #CollisionMesh(eobject.node, name="CollisionMeshAuto1",
                 #surfaceMeshFileName=path+"fingerCollision_part1.stl",
                 #rotation=rotation, translation=translation,
                 #collisionGroup=[1])

    #CollisionMesh(eobject.node, name="CollisionMeshAuto2",
                 #surfaceMeshFileName=path+"fingerCollision_part2.stl",
                 #rotation=rotation, translation=translation,
                 #collisionGroup=[2])
    return finger

class CosseratCable(SofaObject):
    def __init__(self, parentNode, name, trans=[0.0,0.0,0.0], rot=[0.,0.,0.]):
        self.name=name
        self.node = parentNode #.createChild(self.name)             
        self.cableLenght = 81.0
        self.numberBeams = 6
        self.numberFrame = 15
        self.stiffness="50000"
        self.angularStiffness=50000
        self.rot=rot
        self.trans=trans
        self.listOfBeamslenght = None
        self.crossSectionShape='circular'
        self.radius=0.5
        self.youngModulus=5e6    
        
        self.position             = []
        self.frames               = []
        self.cableDofMO           = None
        self.framesDofMO          = []
        self.curv_abs_input       = []
        self.curv_abs_output      = []
        self.cable_position       = []
        self.lenghts              = []
        self.ratePosition         = []
        self.slidingPoint         = None
        self.outputPointMO        = None
        self.mappedPointsNode     = None
        self.framesMO             = None
        
        #self.cable = self.node.createChild(self.name)
        self.__computeFrame()
        self.__computeRate()
        self.__addCables()
    
    def __addCables(self):
        
        rigidBaseNode = self.node.createChild('cableNode')               
        RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", 
                                                 position="0. 0. 0. 0. 0. 0. 1.", translation=self.trans, rotation=self.rot, showObject='1', showObjectScale='0.6',showIndices='1' )
        rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="500", angularStiffness="500",
                                   external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
        
        #############################################
        # Rate of angular Deformation 
        #############################################        
        rateAngularDeformNode = self.node.createChild('rateAngularDeform')
        rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=self.ratePosition)
        BeamHookeLawForce = rateAngularDeformNode.createObject('BeamHookeLawForceField', crossSectionShape=self.crossSectionShape, length=self.lenghts, radius=self.radius, youngModulus=self.youngModulus)
        
        ##############
        #   Frames   #
        ##############
        # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
        mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
        rateAngularDeformNode.addChild(mappedFrameNode)
        framesMO = mappedFrameNode.createObject(
            'MechanicalObject', template='Rigid3d', name="FramesMO", position=self.frames, showObject='1', showObjectScale='1')

        # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
        #                 one output: FramesMO
        inputMO = rateAngularDeformMO.getLinkPath()
        inputMO_rigid = RigidBaseMO.getLinkPath()
        #outputMO = framesMO.getLinkPath()
        self.framesMO = framesMO.getLinkPath()

        curv_abs_input = '0 15 30 45 60 66 81'
        #curv_abs_output = '0.0 5 10 15 20 30 35 40 45 55 60 66 71 76 81'
        mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input=curv_abs_input,
                                    curv_abs_output=self.curv_abs_output, input1=inputMO, input2=inputMO_rigid, output=self.framesMO, debug='0')
        
        #actuators = mappedFrameNode.createChild('actuators')
        #  This create a new node in the scene. This node is appended to the finger's node.
        slidingPoint = mappedFrameNode.createChild('slidingPoint')

        #This create a MechanicalObject, a componant holding the degree of freedom of our
        # mechanical modelling. In the case of a cable it is a set of positions specifying
        #the points where the cable is passing by.
        slidingPointMO = slidingPoint.createObject('MechanicalObject', name="cablePos",
                                position=self.cable_position, showObject="1", showIndices="1")
        slidingPoint.createObject('IdentityMapping')
        
        self.cableDofMO = slidingPointMO.getLinkPath()
        self.slidingPoint = slidingPoint
        
        mappedPointsNode = slidingPoint.createChild('MappedPoints')
        diffPosition = []
        for l in range(0, 7):
            diffPosition.append([0,0,0])
            
        mappedPoints = mappedPointsNode.createObject('MechanicalObject', template='Vec3d', position=diffPosition, name="FramesMO", showObject='1', showObjectScale='1')
        mappedPointsNode.createObject('QPSlidingConstraint', name="QPConstraint")
        
        self.outputPointMO = mappedPoints.getLinkPath()
        self.mappedPointsNode = mappedPointsNode
    
    def __computePosition(self):
        position = []
        for i in range(0,numberBeams):
            position.append([0.,0.,0.])
        
        return position
    
    def __computeRate(self):
        ratePosition = []
        lenghts = []
        
        if self.listOfBeamslenght != None:
            lenghts = self.listOfBeamslenght
        else:
            beamLenght = self.cableLenght/self.numberBeams 
            for i in range(0,self.numberBeams):
                lenghts.append([beamLenght])
                
        for i in range(0,self.numberBeams):
            ratePosition.append([0.,0.,0.])
        
        self.lenghts = lenghts
        self.ratePosition = ratePosition
        
    def __computeFrame(self):
        framLenght = self.cableLenght/self.numberFrame
        
        frame = []
        curv_abs_output = []
        cable_position = []
        for i in range(0,self.numberFrame):
            x = float(i) * framLenght
            curv_abs_output.append(x)
            frame.append([x,0.,0., 0.,0.,0.,1.])
            cable_position.append([x,0.,0.])
            
        frame.append([self.cableLenght,0.,0., 0.,0.,0.,1.])
        curv_abs_output.append(self.cableLenght)
        cable_position.append([self.cableLenght,0.,0.])
        
        self.frames = frame
        self.curv_abs_output = curv_abs_output
        self.cable_position = cable_position
        

    #def constraintBinding(self, position="", attachedTo=None, name = "MappedPoints"):

def CosseratFinger(rootNode,
                   cableNode,
                   translation   =[0., 0., 0.],
                   rotation      =[0., 0., 0.],                   
                   fixingBox        =[-8., -20., -18., 0.0, -3., 8],
                   name             ="1"
                   ):
    cable  = Node(cableNode, name)
    cableN = CosseratCable(cable,
                  name="cable",
                  trans=translation, 
                  rot=rotation)
    
    slidingPoint = cableN.slidingPoint
    cableDofMO   = cableN.cableDofMO
    
    finger = Finger(rootNode, name="Finger"+name, 
                     translation=translation, 
                     rotation=rotation, 
                     fixingBox=fixingBox)
    
    mappedPointsNode = cableN.mappedPointsNode
    inputFEMCableMO = addConstraintPoints(attachedTo=finger,cstPoints= FEMpos,mappedPointsNode=mappedPointsNode, translation=translation,rotation=rotation)
    
    mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, input2=cableDofMO, output=cableN.outputPointMO, direction=cableN.framesMO+".position")
    
    return cable
    
def createScene(rootNode):
    from stlib.scene import MainHeader
    #from stlib.physics.deformable import ElasticMaterialObject

    MainHeader(rootNode, plugins=["CosseratPlugin"], gravity=[0., 0., 0.])
    rootNode.createObject(
        'VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideBoundingCollisionModels hideForceFields showInteractionForceFields showWireframe')
    
    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="5000", printLog="0")
    

    rootNode.gravity = "0 0 0"
    rootNode.createObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.createObject('OglSceneFrame', style="Arrows",alignment="TopRight")
    
    cableNode = rootNode.createChild('cosseratNode')
    cableNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.1')
    cableNode.createObject('SparseLUSolver', name='solver')
    cableNode.createObject('GenericConstraintCorrection')

    cosFinger1 = CosseratFinger(rootNode=rootNode, cableNode=cableNode, 
                         name           ="cosseratF",
                         rotation       =   [-20., 0, 0.],
                         translation    =   [0., 0., 0.0],
                         fixingBox      =   [-20, -20, -10, -2, 10, 15],                         
                         )
    
#    cosseratNode = rootNode.getChild("cosseratNode")
#    cNode            = cosseratNode.getChild("cable")
#    m            = cNode.getChild("cableNode")
    
    GripperController(rootNode, [cosFinger1])
    
#    CF2 = CosseratFinger(rootNode=rootNode, 
#                         cableNode=cableNode, 
#                         name           ="2",
#                         rotation    =[200., 0, 60],
#                         translation =[0., 10., 0.0],
#                         fixingBox      =[-20, -10, 0, 20, 10, 15],
#                         )
#    
#    CF3 = CosseratFinger(rootNode=rootNode, 
#                          cableNode=cableNode, 
#                          name           ="3",
#                          rotation    =[100., 45., 45.],
#                          translation =[-10., 10., -15.0],
#                          fixingBox      =[-20, -10, 0, 20, 10, 15]                          
#                          )
# =============================================================================
    
    #cableNode = rootNode.createChild('cableNode')
    #cableNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.1')
    #cableNode.createObject('SparseLUSolver', name='solver')
    #cableNode.createObject('GenericConstraintCorrection')
    ###### Add cable1 
    #cable1 = CosseratCable(cableNode,name="cable1")
    #slidingPoint = cable1.slidingPoint
    #cableDofMO   = cable1.cableDofMO
    
    #finger1 = Finger(rootNode, translation="-17.5 -12.5 7.5", rotation="0 180 0")
    #mappedPointsNode = cable1.mappedPointsNode
    #inputFEMCableMO = addConstraintPoints(attachedTo=finger1,cstPoints= FEMpos,mappedPointsNode=mappedPointsNode)
    
    #mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, input2=cableDofMO, output=cable1.outputPointMO, direction=cable1.framesMO+".position")
    
    
    
    #cable2 = CosseratCable(cableNode,name="cable2", trans=[0.0,3.0,0.0])
