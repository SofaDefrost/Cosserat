from math import sqrt, pi
from splib.objectmodel import SofaPrefab, SofaObject
from splib.numerics import Vec3, Quat
from splib.animation import animate, AnimationManager
from cosseratUtilities import compute_BeamLenght, createCurvAbsOutput, createFramesList, extractFEMConstraintPoints
import Sofa
import os
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'
dirPath = os.path.dirname(os.path.abspath(__file__))+'/'




def effectorTarget(parentNode, position=[0., 0., 195]):
    target = parentNode.createChild("Target")
    target.createObject("EulerImplicitSolver", firstOrder=True)
    target.createObject("CGLinearSolver")
    target.createObject("MechanicalObject", name="dofs", position=position, showObject=True, showObjectScale=8, drawMode=2, showColor=[1., 1., 1., 1.])
    target.createObject("UncoupledConstraintCorrection")
    return target

## Create all  needs 
## ## ## ## ## ## ## ## 
## Cosserat Needs    ## 
## ## ## ## ## ## ## ## 
length1 = 10.
length2 = 2.
lengthTrunk = 195.
nbBeams = 20
nbFrames = 20 # Dans l'ideal >= nbBeams
nbCables = 4
maxSlidingForce = 0
minSlidingForce = -10

pullPoint = [[0., length1, 0.], [-length1, 0., 0.], [0., -length1, 0.], [length1, 0.,0.]]
direction = Vec3( lengthTrunk, length2-length1, 0.0)
direction.normalize()
Ori0 = Quat.createFromEuler([-pi/2., 0., -0.04], 'ryxz')
Ori1 = Quat.createFromEuler([-pi/2.04, 0.0, 0.0], 'ryxz')
Ori2 = Quat.createFromEuler([-pi/2., 0., 0.04], 'ryxz')
Ori3 = Quat.createFromEuler([-pi/1.96, 0., 0.0], 'ryxz')
baseOrientation = [Ori0,Ori1,Ori2,Ori3]

VecCoord = [] # vector of Vec3
VecFrame = [] # vector of Rigid3d
VecCurvAbsOutput = [] # vector of double
VecCurvAbsInput = [] # vector of double
VecLenght = [] # vector of vetor of double; Cosserat beam length
VecPosFem = []

position = [[0., 0., 0.]]*nbBeams
for k in range(0, nbBeams, 2):
    position[k]   = Vec3(direction[0]*17.5*(k/2)+21, 0.0, 0.0)
    position[k+1] = Vec3(direction[0]*17.5*(k/2)+27, 0.0, 0.0)
position[nbBeams-1][0]  += 1.0 
position=[[0., 0., 0.]]+[pos.toList() for pos in position]
    
####################################################
### Create the Cosserat cable node
####################################################            
VecCoord += position
distance = compute_BeamLenght(position)
VecLenght += distance
VecCurvAbsOutput += createCurvAbsOutput(distance)
VecCurvAbsInput  += createCurvAbsOutput(distance)
VecFrame         += createFramesList(position)
VecPosFem        += extractFEMConstraintPoints(position) # 3D constraint points in the 

#def extractFEMConstraintPoints():    
#    femPoints=[]
#    #### Select point from position
#    for i in range(0,len(position)-1,2):
#        print("=================> i : ",i)
#        femPoints += [[position[i][0], position[i][1], position[i][2]]]                        
#    femPoints += [position[len(position)-1]]
#    
#    #### Apply the predifine transformation    
#    constraintPoints = []
#    for k in range(0, nbCables):
#        constraintPoints.append([])        
#        q = baseOrientation[k]
#        pull = pullPoint[k]
##        for pos in femPoints:
#        for l in range(1, len(femPoints)):
#            pos = femPoints[l]
#            v = Vec3(pos[0], pos[1], pos[2])
#            sol = v.rotateFromQuat(q)
#            #print(" solu :", sol)
#            constraintPoints[k] += [sol[0] + pull[0], sol[1] + pull[1], sol[2] + pull[2]]
#    return constraintPoints

def extractFEMConstraintPoints():    
    femPoints=[]
    #### Select point from position
    for i in range(0,len(position)-1,2):
        center = (position[i][0] + position[i+1][0])/2.0
        femPoints += [[center, position[i][1], position[i][2]]]
    femPoints += [position[len(position)-1]]
    
    #### Apply the predifine transformation    
    constraintPoints = []
    for k in range(0, nbCables):
        constraintPoints.append([])        
        q = baseOrientation[k]
        pull = pullPoint[k]
        for l in range(1, len(femPoints)):
            pos = femPoints[l]
            v = Vec3(pos[0], pos[1], pos[2])
            sol = v.rotateFromQuat(q)
            #print(" solu :", sol)
            constraintPoints[k] += [sol[0] + pull[0], sol[1] + pull[1], sol[2] + pull[2]]
    return constraintPoints

cstPoints = extractFEMConstraintPoints() 

print("Vector of position[0] : ", VecCoord[0])
print("Vector of VecFrame[0] : ", VecFrame[0])


class Animation(Sofa.PythonScriptController):

    def __init__(self, targetNode):
        self.target = targetNode        
        self.rate = 1.
        return

    def initGraph(self, targetNode):
        self.targetMO = self.target.getObject('dofs')
        
    def onKeyPressed(self, c):
        if ord(c) == 21:  # up
            pos = self.targetMO.findData('position').value
            pos[0][0] += self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :", pos)
            
        if ord(c) == 19:  # down
            pos = self.targetMO.findData('position').value
            pos[0][0] -= self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :",pos)
             
        if ord(c) == "+":  # +y
            pos = self.targetMO.findData('position').value
            pos[0][1] += self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :", pos)
            
        if ord(c) == "-":  # -y
            pos = self.targetMO.findData('position').value
            pos[0][1] -= self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :",pos)

        if ord(c) == 18:  # left
            pos = self.targetMO.findData('position').value
            pos[0][2] -= self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :", pos)

        if ord(c) == 20:  # right
            pos = self.targetMO.findData('position').value
            pos[0][2] += self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :", pos)

@SofaPrefab
class CosseratCable(SofaObject):
    def __init__(self, parentNode):
        
        #self.youngModulus = youngModulus
        #self.poissonRatio=poissonRatio        
        self.node = parentNode                
        self.cableDofMOTab = []
        self.outputViolationMOTab = []
        self.mappedPointsNodeTab = []
        self.framesMoTab = []
        
        self.__addCables()
    
    def __addCables(self):
        
        for i in range(0, nbCables):           
            # ###############
            # RigidBases
            ###############
            basePose = VecFrame[0]   
            for l in range(0,3):  ### fill the translation part
                basePose[l] = pullPoint[i][l]
            for l in range(0,4): ### fill the orientation part (here Quat)
                basePose[l+3] = baseOrientation[i][l]
            
            rigidBaseNode = self.node.createChild('rigidBase'+str(i))
            RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", position=basePose, showObject='1', showObjectScale='0.6',showIndices='1' )
            rigidBaseNode.createObject('UniformMass', totalMass="0.001", template="Rigid3d" )
            rigidBaseNode.createObject('PartialFixedConstraint', fixedDirections="1 1 0 1 1 1", indices="0")
#            rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="500", angularStiffness="500", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
            rigidBaseNode.createObject('SlidingActuator', name="SlidingActuator"+str(i), template='Rigid3d', direction='0 0 1 0 0 0', indices=0,  maxForce=maxSlidingForce, minForce=minSlidingForce) 
#            rigidBaseNode.createObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d', direction='1 0 0 0 0 0', indices=0,  maxForce='10000', minForce='-2000') 
            
            #############################################
            # Rate of angular Deformation  (2 sections)
            #############################################
            rate = self.createRateList(position)                           
            longeur = distance  # beams size            
            rateAngularDeformNode = self.node.createChild('rateAngularDeform'+str(i))
            rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO'+str(i), position=rate)
            rateAngularDeformNode.createObject('BeamHookeLawForceField', crossSectionShape='circular', length=longeur, radius='0.2', youngModulus='5e6')
            
            ##############
            #   Frames   #
            ##############
            mappedFrameNode = rigidBaseNode.createChild('MappedFrames'+str(i))
            rateAngularDeformNode.addChild(mappedFrameNode)
            framesMO = mappedFrameNode.createObject('MechanicalObject', template='Rigid3d', name="FramesMO", position=VecFrame, showObject='1', showObjectScale='1', showIndices=1)
            
            inputMO = rateAngularDeformMO.getLinkPath()
            inputMO_rigid = RigidBaseMO.getLinkPath()
            outputMO = framesMO.getLinkPath()            
            mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input=VecCurvAbsInput, curv_abs_output=VecCurvAbsOutput, input1=inputMO, input2=inputMO_rigid, output=outputMO, debug='0', printLog=0) 
            
            ##########################################
            # Multi mapped mstats                    #
            ##########################################
            ###This create a MechanicalObject, a componant holding the degree of freedom of our
            ###mechanical modelling. In the case of a cable it is a set of positions specifying
            ###the points where the cable is passing by.
            slidingPoint = mappedFrameNode.createChild('slidingPoint'+str(i))
            cable_position = self.extractPosfromFrame(framesMO.position)     
            self.framesMoTab += [outputMO]
            slidingPointMO = slidingPoint.createObject('MechanicalObject', name="slidingPointMO"+str(i), position=cable_position, showObject="1", showIndices="1")
            slidingPoint.createObject('IdentityMapping')
            mappedPointsNode = slidingPoint.createChild('MappedPoints')
            mappedPoints = mappedPointsNode.createObject('MechanicalObject', template='Vec3d', position=cstPoints[i], name="FramesMO", showObject='1', showObjectScale='1')
            mappedPointsNode.createObject('CosseratEquality', name="QPConstraint", eqDisp='0.0')
            
            self.mappedPointsNodeTab  += [mappedPointsNode]
            self.cableDofMOTab        += [slidingPointMO.getLinkPath()]                   
            self.outputViolationMOTab += [mappedPoints.getLinkPath()]
   
                        
    def createRateList(self,positions):
        rate=[]
        for i in range(0,len(positions)-1):
            rate += [[0., 0., 0.]]
        return rate
    
    
    
    def getExtractPoint(self):
        print ("END Function getExtractPoint with self.extractPoints =  ",self.extractPoints)
        return self.extractPoints;
    
    def extractPosfromFrame(self,frames):        
        positions = []
        for pos in frames:
            positions += [[pos[0],pos[1],pos[2]]]
        
        return positions
   

@SofaPrefab
class Trunk(SofaObject):
    
    def __init__(self, parentNode, youngModulus=450, poissonRatio=0.45, totalMass=0.042, inverseMode=False):
        
        self.inverseMode = inverseMode
        self.constraintPointMoTab = []
        self.node = parentNode.createChild('Trunk')

        self.node.createObject('MeshVTKLoader', name='loader', filename=path+'trunk.vtk')
        self.node.createObject('TetrahedronSetTopologyContainer', src='@loader', name='container')
        self.node.createObject('TetrahedronSetTopologyModifier')
        self.node.createObject('TetrahedronSetTopologyAlgorithms')
        self.node.createObject('TetrahedronSetGeometryAlgorithms')

        self.node.createObject('MechanicalObject', name='dofs', template='Vec3d', showIndices='false', showIndicesScale='4e-5')
        self.node.createObject('UniformMass', totalMass=totalMass)
        self.node.createObject('TetrahedronFEMForceField', template='Vec3d', name='FEM', method='large', poissonRatio=poissonRatio,  youngModulus=youngModulus)

        
    def addVisualModel(self, color=[1., 1., 1., 1.]):
        trunkVisu = self.node.createChild('VisualModel')
        trunkVisu.createObject('MeshSTLLoader', filename=path+"trunk.stl", rotation="0 90 0")
        trunkVisu.createObject('OglModel', template='ExtVec3d', color=color)
        trunkVisu.createObject('BarycentricMapping')
        

    def addCollisionModel(self, selfCollision=False):
        trunkColli = self.node.createChild('CollisionModel')
        for i in range(2):
            part = trunkColli.createChild("Part"+str(i+1))
            part.createObject('MeshSTLLoader', name="loader", filename=path+"trunk_colli"+str(i+1)+".stl")
            part.createObject('MeshTopology', src="@loader")
            part.createObject('MechanicalObject', rotation="0 90 0")
            part.createObject('TTriangleModel', group=1 if not selfCollision else i)
            part.createObject('TLineModel', group=1 if not selfCollision else i)
            part.createObject('TPointModel', group=1 if not selfCollision else i)
            part.createObject('BarycentricMapping')

    def fixExtremity(self):
        self.node.createObject('BoxROI', name='boxROI', box=[[-20, -20, 0], [20, 20, 20]], drawBoxes=False)
        self.node.createObject('PartialFixedConstraint', fixedDirections="1 1 1", indices="@boxROI.indices")

    def addEffectors(self, target, position=[0., 0., 195.]):
        effectors = self.node.createChild("Effectors")
        effectors.createObject("MechanicalObject", position=position)
        effectors.createObject("PositionEffector", indices=range(len(position)), effectorGoal=target)
        effectors.createObject("BarycentricMapping", mapForces=False, mapMasses=False)
    
    ##########################################
    #       add constraint points            #
    ##########################################
    def addConstraintPoints(self, cstPoints,i,mappedPointsNode):
        trunkMappedPoints = self.node.createChild('trunkMappedPoints'+str(i))        
        inputFEMCable = trunkMappedPoints.createObject('MechanicalObject', name="pointsInFEM"+str(i), position=cstPoints, showObject="1", showIndices="1")
        trunkMappedPoints.addChild(mappedPointsNode)
        trunkMappedPoints.createObject('BarycentricMapping')
        
        self.constraintPointMoTab += [inputFEMCable.getLinkPath()]


def createScene(rootNode):

    rootNode.createObject("RequiredPlugin", name="SoftRobots")
    rootNode.createObject("RequiredPlugin", name="SoftRobots.Inverse")
    rootNode.createObject("RequiredPlugin", name="SofaSparseSolver")
    rootNode.createObject("RequiredPlugin", name="SofaPreconditioner")
    rootNode.createObject("RequiredPlugin", name="SofaPython")
    rootNode.createObject("RequiredPlugin", name="CosseratPlugin")
    
    AnimationManager(rootNode)
    rootNode.createObject("VisualStyle", displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideBoundingCollisionModels hideForceFields showInteractionForceFields showWireframe')
    rootNode.gravity = "0 0 0"

    rootNode.createObject("FreeMotionAnimationLoop")
    #### For direct resolution, i.e direct control of the cable displacement    
    # rootNode.createObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")
    #### For inverse resolution, i.e control of effectors position
    rootNode.createObject("QPInverseProblemSolver", printLog='0', epsilon=1e-1, maxIterations="500")

    # ###############
    # New adds to use the sliding Actuator
    ###############
    cableNode = rootNode.createChild('cableNode')
    cableNode.createObject('EulerImplicitSolver', firstOrder="0",rayleighMass="0.1", rayleighStiffness="0.1")
    cableNode.createObject('SparseLUSolver', name='solver')
    #cableNode.createObject('SparseLDLSolver', name='solver')
    cableNode.createObject('GenericConstraintCorrection')
    
    #tabOfNode = []
    Cable = CosseratCable(cableNode)       
    simulation = rootNode.createChild("Simulation")

    simulation.createObject('EulerImplicitSolver', name='odesolver', firstOrder="0", rayleighMass="0.1", rayleighStiffness="0.1")
    simulation.createObject('ShewchukPCGLinearSolver', name='linearSolver', iterations='500', tolerance='1.0e-18', preconditioners="precond")
    simulation.createObject('SparseLDLSolver', name='precond')
    simulation.createObject('GenericConstraintCorrection', solverName="precond")

    trunk = Trunk(simulation, inverseMode=True)
    trunk.addVisualModel(color=[1., 1., 1., 0.8])
    trunk.fixExtremity()
    
    
    mappedPointsNodeTab = Cable.mappedPointsNodeTab
    cableDofMOTab = Cable.cableDofMOTab
    outputViolationMOTab = Cable.outputViolationMOTab
    framesMoTab = Cable.framesMoTab
    for i in range(0,nbCables):
        mappedPointsNode = mappedPointsNodeTab[i]
        trunk.addConstraintPoints(cstPoints[i],i,mappedPointsNode)
        
        #### Get link to different Mo for the multi map
        constraintPointMo = trunk.constraintPointMoTab[i]
        cableDofMO = cableDofMOTab[i]
        outputViolationMO = outputViolationMOTab[i]
        print(" ==================+++++> The Link is : ", framesMoTab[i])
        
        mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=constraintPointMo, input2=cableDofMO, output=outputViolationMO, direction=framesMoTab[i]+".position")
        
    
    target = effectorTarget(rootNode)
    trunk.addEffectors(target=target.dofs.getData("position").getLinkPath(), position=[[0., 0., 195]])


    ################################
    # Animation (to move the dofs) #
    ################################
    Animation(target)
    
    #### Use this in direct mode as an example of animation ############
    #    def cableanimation(target, factor):
    #        target.cable.value = factor*20
    
    #    animate(cableanimation, {"target": trunk.cableL0}, duration=2, )
    #################################################################
