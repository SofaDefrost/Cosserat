# coding: utf8

import Sofa
import Sofa.Simulation
import SofaRuntime
import unittest
import numpy

class Test(unittest.TestCase):
    def test_example_scene(self):

        scene = Sofa.import_sofa_python_scene("/home/dmarchal/projects/dev/sofa1/plugins/Cosserat/docs/testScene/tuto_1.py")

        # create a root node to store the scene
        root = Sofa.Core.Node("rootNode")

        # fill the root node with the content of the scene
        scene.createScene(root)

        # initialize everything
        Sofa.Simulation.init(root)

        scenarios = [{
                "value" : [[0.0,0.0,0.0]] * 3,
                "result" :  [[-0., -0., -0.,  0.,  0.,  0.,  1.],
                             [10., -0., -0.,  0.,  0.,  0.,  1.],
                             [20., -0., -0.,  0.,  0.,  0.,  1.],
                             [30., -0., -0.,  0.,  0.,  0.,  1.]]
            },
            {
                "value" : [[0.0,0.1,0.0]] * 3,
                "result" :  [[ -0.        ,  -0.        ,  -0.        ,   0.        ,
                                0.        ,   0.        ,   1.        ],
                                [  8.41470985,  -0.        ,  -4.59697694,   0.        ,
                                   0.47942554,   0.        ,   0.87758256],
                                [  9.09297427,  -0.        , -14.16146837,   0.        ,
                                   0.84147098,   0.        ,   0.54030231],
                                [  1.41120008,  -0.        , -19.89992497,   0.        ,
                                    0.99749499,   0.        ,   0.0707372 ]]
            }
        ]
        for scenario in scenarios:
            root.cosseratCoordinate.cosserat_state.position.value = scenario["value"]
            Sofa.Simulation.animate(root, 0.01)
            numpy.testing.assert_array_almost_equal(root.rigid_base.cosserat_in_Sofa_frame_node.FramesMO.position.value, scenario["result"])

