from typing import Any, List, Optional, Tuple, Union

import numpy as np
import Sofa


def addEdgeCollision(parentNode: Sofa.Core.Node,
                     position3D: List[List[float]],
                     edges: List[List[int]],
                     group: str = '2') -> Sofa.Core.Node:
    """
    Add edge-based collision model to a parent node.

    This function creates a child node with edge collision models for detecting
    collisions between a Cosserat rod and other objects in the scene.

    Args:
        parentNode: The parent node to attach the collision model to
        position3D: List of 3D positions for collision vertices
        edges: List of edge indices connecting vertices
        group: Collision group identifier (default: '2')

    Returns:
        The created collision node

    Example:
        ```python
        # Create collision model
        edges = [[0, 1], [1, 2], [2, 3]]
        positions = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]
        collision_node = addEdgeCollision(beam_node, positions, edges)
        ```
    """
    if not parentNode:
        raise ValueError("Parent node cannot be None")

    if not position3D:
        raise ValueError("position3D cannot be empty")

    if not edges:
        raise ValueError("edges cannot be empty")

    collisInstrumentCombined = parentNode.addChild('collisInstrumentCombined')
    collisInstrumentCombined.addObject('EdgeSetTopologyContainer',
                                      name="collisEdgeSet",
                                      position=position3D,
                                      edges=edges)
    collisInstrumentCombined.addObject('EdgeSetTopologyModifier',
                                      name="collisEdgeModifier")
    collisInstrumentCombined.addObject('MechanicalObject',
                                      name="CollisionDOFs")
    collisInstrumentCombined.addObject('LineCollisionModel',
                                      bothSide="1",
                                      group=group)
    collisInstrumentCombined.addObject('PointCollisionModel',
                                      group=group)
    collisInstrumentCombined.addObject('IdentityMapping',
                                      name="mapping")
    return collisInstrumentCombined


def addPointsCollision(parentNode: Sofa.Core.Node,
                      position3D: List[List[float]],
                      edges: List[List[int]],
                      nodeName: str,
                      group: str = '2') -> Sofa.Core.Node:
    """
    Add point-based collision model to a parent node.

    This function creates a child node with point collision models, which is
    especially useful for detecting interactions at beam endpoints or specific points.

    Args:
        parentNode: The parent node to attach the collision model to
        position3D: List of 3D positions for collision vertices
        edges: List of edge indices connecting vertices
        nodeName: Name for the created collision node
        group: Collision group identifier (default: '2')

    Returns:
        The created collision node

    Example:
        ```python
        # Create point collision model
        edges = [[0, 1], [1, 2], [2, 3]]
        positions = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]
        collision_node = addPointsCollision(beam_node, positions, edges, "BeamCollision")
        ```
    """
    if not parentNode:
        raise ValueError("Parent node cannot be None")

    if not position3D:
        raise ValueError("position3D cannot be empty")

    if not edges:
        raise ValueError("edges cannot be empty")

    if not nodeName:
        raise ValueError("nodeName cannot be empty")

    collisInstrumentCombined = parentNode.addChild(nodeName)
    collisInstrumentCombined.addObject('EdgeSetTopologyContainer',
                                      name="beamContainer",
                                      position=position3D,
                                      edges=edges)
    collisInstrumentCombined.addObject('EdgeSetTopologyModifier',
                                      name="beamModifier")
    collisInstrumentCombined.addObject('MechanicalObject',
                                      name="collisionStats",
                                      showObject=False,
                                      showIndices=False)
    collisInstrumentCombined.addObject('PointCollisionModel',
                                      name="beamColMod",
                                      group=group)
    # Using RigidMapping instead of IdentityMapping for better performance with rigid bodies
    collisInstrumentCombined.addObject('RigidMapping',
                                      name="beamMapping")
    return collisInstrumentCombined


def addConstraintPoint(parentNode: Sofa.Core.Node,
                    beamPath: str = "/solverNode/needle/rigidBase/cosseratInSofaFrameNode/slidingPoint/slidingPointMO") -> Sofa.Core.Node:
    """
    Build a constraint node for applying constraints to a Cosserat rod.

    This function creates a constraint points node that can be used to apply
    constraints at specific points along a beam.

    Args:
        parentNode: The parent node to attach the constraint model to
        beamPath: Path to the beam's mechanical object (default points to a standard needle path)

    Returns:
        The created constraint node

    Example:
        ```python
        # Create constraint node
        beam_path = "/path/to/beam/mechanicalObject"
        constraint_node = addConstraintPoint(root_node, beam_path)
        ```
    """
    if not parentNode:
        raise ValueError("Parent node cannot be None")

    constraintPointsNode = parentNode.addChild('constraintPoints')
    constraintPointsNode.addObject("PointSetTopologyContainer",
                                  name="constraintPtsContainer",
                                  listening="1")
    constraintPointsNode.addObject("PointSetTopologyModifier",
                                  name="constraintPtsModifier",
                                  listening="1")
    constraintPointsNode.addObject("MechanicalObject",
                                  template="Vec3d",
                                  showObject=True,
                                  showIndices=True,
                                  name="constraintPointsMo",
                                  position=[],
                                  showObjectScale=0,
                                  listening="1")

    constraintPointsNode.addObject('PointsManager',
                                  name="pointsManager",
                                  listening="1",
                                  beamPath=beamPath)

    constraintPointsNode.addObject('BarycentricMapping',
                                  useRestPosition="false",
                                  listening="1")
    return constraintPointsNode


def addSlidingPoints(parentNode: Sofa.Core.Node,
                    frames3D: List[List[float]],
                    showVisual: bool = False) -> Sofa.Core.Node:
    """
    Add sliding points to a parent node for sliding contact simulation.

    This function creates a child node with sliding points that can be used
    to represent points that slide along another object.

    Args:
        parentNode: The parent node to attach the sliding points to
        frames3D: List of 3D positions for the sliding points
        showVisual: Whether to show the points visually (default: False)

    Returns:
        The created sliding points node

    Example:
        ```python
        # Create sliding points
        positions = [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]
        sliding_node = addSlidingPoints(beam_node, positions)
        ```
    """
    if not parentNode:
        raise ValueError("Parent node cannot be None")

    if not frames3D:
        raise ValueError("frames3D cannot be empty")

    slidingPoint = parentNode.addChild('slidingPoint')
    slidingPoint.addObject('MechanicalObject',
                          name="slidingPointMO",
                          position=frames3D,
                          showObject=showVisual,
                          showIndices=showVisual)
    slidingPoint.addObject('IdentityMapping')
    return slidingPoint


def getLastConstraintPoint(constraintPointsNode: Sofa.Core.Node) -> np.ndarray:
    """
    Get the last constraint point position from a constraint node.

    Args:
        constraintPointsNode: The constraint points node

    Returns:
        The position of the last constraint point as a numpy array [x, y, z]

    Raises:
        ValueError: If the node doesn't have a 'constraintPointsMo' object or if there are no points

    Example:
        ```python
        last_point = getLastConstraintPoint(constraint_node)
        print(f"Last constraint point is at: {last_point}")
        ```
    """
    if not constraintPointsNode:
        raise ValueError("Constraint points node cannot be None")

    try:
        mo = constraintPointsNode.getObject('constraintPointsMo')
        if len(mo.position) == 0:
            raise ValueError("No constraint points available")
        return mo.position[-1]
    except Exception as e:
        raise ValueError(f"Error accessing constraint points: {e}")


def computeDistanceBetweenPoints(constraintPointPos: List[np.ndarray],
                               slidingPointPos: List[np.ndarray]) -> float:
    """
    Compute the Euclidean distance between the last constraint point and sliding point.

    This function calculates the 3D distance between the last points in the provided arrays,
    which is useful for determining proximity or contact between points.

    Args:
        constraintPointPos: List of constraint point positions as numpy arrays
        slidingPointPos: List of sliding point positions as numpy arrays

    Returns:
        The distance between the last points, or 0 if no constraint points exist

    Example:
        ```python
        constraint_points = [[0, 0, 0], [1, 0, 0]]
        sliding_points = [[0, 1, 0], [1, 1, 0]]
        distance = computeDistanceBetweenPoints(constraint_points, sliding_points)
        print(f"Distance: {distance}")  # Output: Distance: 1.0
        ```
    """
    if len(constraintPointPos) == 0:
        return 0.0

    try:
        return np.linalg.norm(constraintPointPos[-1] - slidingPointPos[-1])
    except (IndexError, ValueError) as e:
        print(f"Error computing distance: {e}")
        return 0.0


def computePositiveAlongXDistanceBetweenPoints(constraintPointPos: List[np.ndarray],
                                            slidingPointPos: List[np.ndarray]) -> float:
    """
    Compute the distance between points only if the constraint point is ahead along X-axis.

    This function calculates the distance only when the constraint point has a greater
    X-coordinate than the sliding point, otherwise returns 0.

    Args:
        constraintPointPos: List of constraint point positions as numpy arrays
        slidingPointPos: List of sliding point positions as numpy arrays

    Returns:
        The distance between the last points if constraint is ahead in X, otherwise 0
    """
    if len(constraintPointPos) == 0:
        return 0.0

    try:
        if constraintPointPos[-1][0] > slidingPointPos[-1][0]:
            return np.linalg.norm(constraintPointPos[-1] - slidingPointPos[-1])
        else:
            return 0.0
    except (IndexError, ValueError) as e:
        print(f"Error computing positive X distance: {e}")
        return 0.0


def computeNegativeAlongXDistanceBetweenPoints(constraintPointPos: List[np.ndarray],
                                            slidingPointPos: List[np.ndarray]) -> float:
    """
    Compute the distance between points only if the constraint point is behind along X-axis.

    This function calculates the distance only when the constraint point has a smaller
    X-coordinate than the sliding point, otherwise returns 0.

    Args:
        constraintPointPos: List of constraint point positions as numpy arrays
        slidingPointPos: List of sliding point positions as numpy arrays

    Returns:
        The distance between the last points if constraint is behind in X, otherwise 0
    """
    if len(constraintPointPos) == 0:
        return 0.0

    try:
        if constraintPointPos[-1][0] < slidingPointPos[-1][0]:
            return np.linalg.norm(constraintPointPos[-1] - slidingPointPos[-1])
        else:
            return 0.0
    except (IndexError, ValueError) as e:
        print(f"Error computing negative X distance: {e}")
        return 0.0


def create_rigid_node(parent_node: Sofa.Core.Node,
                     name: str,
                     translation = [0.0, 0.0, 0.0] ,
                     rotation = [0.0, 0.0, 0.0] ,
                     positions: List[List[float]] = None) -> Sofa.Core.Node:
    """
    Create a rigid body node with mechanical object.

    This function creates a child node with a rigid body mechanical object,
    which is useful for representing rigid parts of a Cosserat rod or for
    attaching other objects to the rod.

    Args:
        parent_node: The parent node to attach the rigid node to
        name: Name for the created rigid node
        translation: Initial translation [x, y, z]
        rotation: Initial rotation [rx, ry, rz] (Euler angles in radians)
        positions: Initial positions of the rigid body [tx, ty, tz, rx, ry, rz, w]
                  (default: [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])

    Returns:
        The created rigid node

    Example:
        ```python
        # Create rigid base for a rod
        translation = [0, 0, 0]
        rotation = [0, 0, 0]
        rigid_node = create_rigid_node(root_node, "RigidBase", translation, rotation)
        ```
    """
    if not parent_node:
        raise ValueError("Parent node cannot be None")

    if not name:
        raise ValueError("Name cannot be empty")

    if positions is None:
        positions = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]

    # Validate translation and rotation
    if not isinstance(translation, (list, tuple, np.ndarray)) or len(translation) != 3:
        raise ValueError("Translation must be a list of 3 values [x, y, z]")

    if not isinstance(rotation, (list, tuple, np.ndarray)) or len(rotation) != 3:
        raise ValueError("Rotation must be a list of 3 values [rx, ry, rz]")

    rigidBaseNode = parent_node.addChild(name)
    rigidBaseNode.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name=f"{name}MO",
        position=positions,
        translation=translation,
        rotation=rotation
    )

    return rigidBaseNode


