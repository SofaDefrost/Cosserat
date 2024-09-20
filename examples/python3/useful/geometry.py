#
from typing import List
from useful.params import BeamGeometryParameters

def calculate_beam_parameters(beamGeoParams):
    # Data validation checks for beamGeoParams attributes
    if not all(hasattr(beamGeoParams, attr) for attr in ['init_pos', 'beamLength', 'nbSection']):
        raise ValueError("beamGeoParams must have 'init_pos', 'beamLength', and 'nbSection' attributes.")

    total_length = beamGeoParams.beamLength
    nb_sections = beamGeoParams.nbSection
    x, y, z = beamGeoParams.init_pos

    if not all(isinstance(val, (int, float)) for val in [x, y, z, total_length]):
        raise ValueError("init_pos and beamLength in beamGeoParams must be numeric values.")

    if not isinstance(nb_sections, int) or nb_sections <= 0:
        raise ValueError("nbSection in beamGeoParams must be a positive integer.")

    length_s = total_length / nb_sections
    bendingState = []
    listOfSectionsLength = []
    temp = x
    curv_abs_input_s = [x]

    for i in range(nb_sections):
        bendingState.append([0, 0, 0])
        listOfSectionsLength.append((((i + 1) * length_s) - i * length_s))
        temp += listOfSectionsLength[i]
        curv_abs_input_s.append(temp)
    curv_abs_input_s[nb_sections] = total_length + x

    return bendingState, curv_abs_input_s, listOfSectionsLength


def calculate_frame_parameters(beamGeoParams):
    # Data validation checks for beamGeoParams attributes
    if not all(hasattr(beamGeoParams, attr) for attr in ['init_pos', 'beamLength', 'nbFrames']):
        raise ValueError("beamGeoParams must have 'init_pos', 'beamLength', and 'nbFrames' attributes.")

    x, y, z = beamGeoParams.init_pos
    total_length = beamGeoParams.beamLength
    nb_frames = beamGeoParams.nbFrames

    if not all(isinstance(val, (int, float)) for val in [x, y, z, total_length]):
        raise ValueError("init_pos and beamLength in beamGeoParams must be numeric values.")

    if not isinstance(nb_frames, int) or nb_frames <= 0:
        raise ValueError("nbFrames in beamGeoParams must be a positive integer.")

    length_f = total_length / nb_frames
    frames_f = []
    curv_abs_output_f = []
    cable_position_f = []

    # @Todo: improve this for
    for i in range(nb_frames):
        sol = i * length_f
        frames_f.append([sol + x, y, z, 0, 0, 0, 1])
        cable_position_f.append([sol + x, y, z])
        curv_abs_output_f.append(sol + x)

    frames_f.append([total_length + x, y, z, 0, 0, 0, 1])
    cable_position_f.append([total_length + x, y, z])
    curv_abs_output_f.append(total_length + x)

    return frames_f, curv_abs_output_f, cable_position_f


def generate_edge_list(cable3DPos: List[List[float]]) -> list[list[int]]:
    """
    Generate an edge list required in the EdgeSetTopologyContainer component.

    Parameters:
        cable3DPos (List[List[float]]): A list of 3D points representing the cable positions.

    Returns:
        List[int]: A list of indices forming edges in the EdgeSetTopologyContainer.
    """
    number_of_points = len(cable3DPos)
    edges = []
    for i in range(number_of_points - 1):
        edges.append([i,i+1])
    return edges



class CosseratGeometry:
    def __init__(self, beamGeoParams):
        # Data validation checks for beamGeoParams
        if not isinstance(beamGeoParams, BeamGeometryParameters):
            raise ValueError("beamGeoParams must be an instance of BeamGeoParams.")

        self.bendingState, self.curv_abs_inputS, self.sectionsLengthList = calculate_beam_parameters(beamGeoParams)
        self.framesF, self.curv_abs_outputF, self.cable_positionF = calculate_frame_parameters(beamGeoParams)
