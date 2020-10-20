from math import cos
from math import sin, sqrt, pi
from splib.objectmodel import SofaPrefab, SofaObject
from splib.numerics import Vec3, Quat
from splib.animation import animate, AnimationManager

import os
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'
dirPath = os.path.dirname(os.path.abspath(__file__))+'/'


def effectorTarget(parentNode, position=[0., 0., 200]):
    target = parentNode.createChild("Target")
    target.createObject("EulerImplicitSolver", firstOrder=True)
    target.createObject("CGLinearSolver")
    target.createObject("MechanicalObject", name="dofs", position=position, showObject=True, showObjectScale=8, drawMode=2, showColor=[1., 1., 1., 1.])
    target.createObject("UncoupledConstraintCorrection")
    return target


@SofaPrefab
class Trunk(SofaObject):
    """ This prefab is implementing a soft robot inspired by the elephant's trunk.
        The robot is entirely soft and actuated with 8 cables.

        The prefab is composed of:
        - a visual model
        - a collision model
        - a mechanical model for the deformable structure

        The prefab has the following parameters:
        - youngModulus
        - poissonRatio
        - totalMass

        Example of use in a Sofa scene:

        def createScene(root):
            ...
            trunk = Trunk(root)

            ## Direct access to the components
            trunk.displacements = [0., 0., 0., 0., 5., 0., 0., 0.]
    """

    def __init__(self, parentNode, youngModulus=450, poissonRatio=0.45, totalMass=0.042, inverseMode=False):

        self.inverseMode = inverseMode
        self.node = parentNode.createChild('Trunk')

        self.node.createObject('MeshVTKLoader', name='loader', filename=path+'trunk.vtk')
        self.node.createObject('TetrahedronSetTopologyContainer', src='@loader', name='container')
        self.node.createObject('TetrahedronSetTopologyModifier')
        self.node.createObject('TetrahedronSetTopologyAlgorithms')
        self.node.createObject('TetrahedronSetGeometryAlgorithms')

        self.node.createObject('MechanicalObject', name='dofs', template='Vec3d', showIndices='false', showIndicesScale='4e-5')
        self.node.createObject('UniformMass', totalMass=totalMass)
        self.node.createObject('TetrahedronFEMForceField', template='Vec3d', name='FEM', method='large', poissonRatio=poissonRatio,  youngModulus=youngModulus)

        self.__addCables()
        
    def compute_distance(self, positions):
        listDist=[]
        listDist += [0]
        for i in range(0,len(positions)-1):
            t1 = positions[i]
            t2 = positions[i+1]
            dist = sqrt((t1[0]-t2[0])**2 + (t1[1]-t2[1])**2 + (t1[2]-t2[2])**2)
            
            listDist += [dist] 
            
        return listDist
    
    def createFramesList(self,positions):
        frames=[]
        for i in range(0,len(positions)):
            frames += [[positions[i][0],positions[i][1],positions[i][2],0,0,0,1]]
            
        return frames
    
    def createCurvAbsOutput(self,distance):
        curv_abs_output=[]
        
        for i in range(0,len(distance)):
            temp = 0
            for k in range(0,i+1):
                temp += distance[k]
            curv_abs_output += [temp]
        
        return curv_abs_output
    
    
    
    def __addCables(self):
        length1 = 10.
        length2 = 2.
        lengthTrunk = 195.

        pullPoint = [[0., length1, 0.], [-length1, 0., 0.], [0., -length1, 0.], [length1, 0.,0.]]
        direction = Vec3( lengthTrunk, length2-length1, 0.0)
        direction.normalize()
        Ori0 = Quat.createFromEuler([-pi/2., 0., -0.04], 'ryxz')
        Ori1 = Quat.createFromEuler([-pi/2.04, 0.0, 0.0], 'ryxz')
        Ori2 = Quat.createFromEuler([-pi/2., 0., 0.04], 'ryxz')
        Ori3 = Quat.createFromEuler([-pi/1.96, 0., 0.0], 'ryxz')
        
        baseOrientation = [Ori0,Ori1,Ori2,Ori3]
        
        nbCables = 4
        for i in range(0, nbCables):
            position = [[0., 0., 0.]]*20
            for k in range(0, 20, 2):
                position[k]   = Vec3(direction[0]*17.5*(k/2)+21, 0.0, 0.0)                 
                position[k+1] = Vec3(direction[0]*17.5*(k/2)+27, 0.0, 0.0)
            
            position=[[0., 0, 0.]]+[pos.toList() for pos in position]
            
            ####################################################
            ### Create the Cosserat cable node
            ####################################################            
            distance = self.compute_distance(position)            
            curv_abs_output = self.createCurvAbsOutput(distance)
            curv_abs_input = self.createCurvAbsOutput(distance)
            frames = self.createFramesList(position)
            
            # ###############
            # RigidBases
            ###############
            basePose = frames[0]            
            for l in range(0,3):
                basePose[l] = pullPoint[i][l]
            for l in range(0,4):
                basePose[l+3] = baseOrientation[i][l]
            
            rigidBaseNode = self.node.createChild('rigidBase'+str(i))
            RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", position=basePose, showObject='1', showObjectScale='0.6',showIndices='1' )
            rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="500", angularStiffness="500", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
            
            #############################################
            # Rate of angular Deformation  (2 sections)
            #############################################
            ratPosition = []
            for l in range(0,len(frames)-1):
                ratPosition += [0,0,0]
            longeur = distance  # beams size
            rateAngularDeformNode = self.node.createChild('rateAngularDeform'+str(i))
            rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=ratPosition)
            BeamHookeLawForce = rateAngularDeformNode.createObject('BeamHookeLawForceField', crossSectionShape='circular', length=longeur, radius='0.5', youngModulus='5e6',showIndices='1')
            
            ##############
            #   Frames   #
            ##############
            mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
            rateAngularDeformNode.addChild(mappedFrameNode)
            framesMO = mappedFrameNode.createObject('MechanicalObject', template='Rigid3d', name="FramesMO", position=frames, showObject='1', showObjectScale='1', showIndices=1)
            
            inputMO = rateAngularDeformMO.getLinkPath()
            inputMO_rigid = RigidBaseMO.getLinkPath()
            outputMO = framesMO.getLinkPath()
            
            mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input=curv_abs_input, curv_abs_output=curv_abs_output, input1=inputMO, input2=inputMO_rigid, output=outputMO, debug='0', printLog=0) 
        
    def addVisualModel(self, color=[1., 1., 1., 1.]):
        trunkVisu = self.node.createChild('VisualModel')
        trunkVisu.createObject('MeshSTLLoader', filename=path+"trunk.stl", rotation="0 0 0")
        trunkVisu.createObject('OglModel', template='ExtVec3d', color=color)
        trunkVisu.createObject('BarycentricMapping')

    def addCollisionModel(self, selfCollision=False):
        trunkColli = self.node.createChild('CollisionModel')
        for i in range(2):
            part = trunkColli.createChild("Part"+str(i+1))
            part.createObject('MeshSTLLoader', name="loader", filename=path+"trunk_colli"+str(i+1)+".stl")
            part.createObject('MeshTopology', src="@loader")
            part.createObject('MechanicalObject', rotation="0 0 0")
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


def createScene(rootNode):

    rootNode.createObject("RequiredPlugin", name="SoftRobots")
    rootNode.createObject("RequiredPlugin", name="SoftRobots.Inverse")
    rootNode.createObject("RequiredPlugin", name="SofaSparseSolver")
    rootNode.createObject("RequiredPlugin", name="SofaPreconditioner")
    rootNode.createObject("RequiredPlugin", name="SofaPython")
    rootNode.createObject("RequiredPlugin", name="CosseratPlugin")
    
    AnimationManager(rootNode)
    rootNode.createObject("VisualStyle", displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideBoundingCollisionModels hideForceFields showInteractionForceFields showWireframe')
    rootNode.gravity = [0., -9810., 0.]

    rootNode.createObject("FreeMotionAnimationLoop")
    # For direct resolution, i.e direct control of the cable displacement
    # rootNode.createObject("GenericConstraintSolver", maxIterations=100, tolerance=1e-5)
    # For inverse resolution, i.e control of effectors position
    rootNode.createObject("QPInverseProblemSolver", epsilon=1e-1)

    simulation = rootNode.createChild("Simulation")

    simulation.createObject('EulerImplicitSolver', name='odesolver', firstOrder="0", rayleighMass="0.1", rayleighStiffness="0.1")
    simulation.createObject('ShewchukPCGLinearSolver', name='linearSolver', iterations='500', tolerance='1.0e-18', preconditioners="precond")
    simulation.createObject('SparseLDLSolver', name='precond')
    simulation.createObject('GenericConstraintCorrection', solverName="precond")

    trunk = Trunk(simulation, inverseMode=True)
    trunk.addVisualModel(color=[1., 1., 1., 0.8])
    trunk.fixExtremity()
    target = effectorTarget(rootNode)
    trunk.addEffectors(target=target.dofs.getData("position").getLinkPath(), position=[[0., 0., 195]])

    # Use this in direct mode as an example of animation ############
    # def cableanimation(target, factor):
    #     target.cable.value = factor*20
    #
    # animate(cableanimation, {"target": trunk.cableL0}, duration=2, )
    #################################################################
