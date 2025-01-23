#
from useful.params import BeamGeometryParameters


def calculate_beam_parameters(beamGeoParams):
    # Data validation checks for beamGeoParams attributes
    if not all(hasattr(beamGeoParams, attr) for attr in ['beam_length', 'nb_section']):
        raise ValueError("beamGeoParams must have, 'beam_length', and 'nb_section' attributes.")

    total_length = beamGeoParams.beam_length
    nb_sections = beamGeoParams.nb_section

    if not all(isinstance(val, (float)) for val in [total_length]):
        raise ValueError("init_pos and beamLength in beamGeoParams must be numeric values.")

    if not isinstance(nb_sections, int) or nb_sections <= 0:
        raise ValueError("nbSection in beamGeoParams must be a positive integer.")

    length_s = total_length / nb_sections
    bendingState = []
    listOfSectionsLength = []
    temp = 0
    curv_abs_input_s = [0]

    for i in range(nb_sections):
        bendingState.append([0, 0, 0])
        listOfSectionsLength.append((((i + 1) * length_s) - i * length_s))
        temp += listOfSectionsLength[i]
        curv_abs_input_s.append(temp)
    curv_abs_input_s[nb_sections] = total_length

    return bendingState, curv_abs_input_s, listOfSectionsLength


def calculate_frame_parameters(beamGeoParams):
    # Data validation checks for beamGeoParams attributes
    if not all(hasattr(beamGeoParams, attr) for attr in ['beam_length', 'nb_frames']):
        raise ValueError("beamGeoParams must have 'beamLength', and 'nbFrames' attributes.")

    total_length = beamGeoParams.beam_length
    nb_frames = beamGeoParams.nb_frames

    if not all(isinstance(val, (int, float)) for val in [total_length]):
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
        frames_f.append([sol, 0, 0, 0, 0, 0, 1])
        cable_position_f.append([sol, 0, 0])
        curv_abs_output_f.append(sol)

    frames_f.append([total_length, 0, 0, 0, 0, 0, 1])
    cable_position_f.append([total_length, 0, 0])
    curv_abs_output_f.append(total_length)

    return frames_f, curv_abs_output_f, cable_position_f


def generate_edge_list(cable3DPos: list[list[float]]) -> list[list[int]]:
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
        edges.append([i, i + 1])
    return edges


class CosseratGeometry:
    def __init__(self, beamGeoParams):
        # Data validation checks for beamGeoParams
        if not isinstance(beamGeoParams, BeamGeometryParameters):
            raise ValueError("beamGeoParams must be an instance of BeamGeoParams.")

        self.bendingState, self.curv_abs_inputS, self.sectionsLengthList = calculate_beam_parameters(beamGeoParams)
        self.framesF, self.curv_abs_outputF, self.cable_positionF = calculate_frame_parameters(beamGeoParams)
