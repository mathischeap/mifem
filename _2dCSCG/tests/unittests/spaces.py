
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from _2dCSCG.master import SpaceInvoker



def test_Space_NO1_polynomial_space():
    """ Unittests for the mesh."""
    if rAnk == mAster_rank:
        print(">>> [test_Space_NO1_polynomial_space] ...... ", flush=True)

    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 2)])

    np.testing.assert_array_almost_equal(space.nodes[0],
                                         np.array([-1.       , -0.4472136,  0.4472136,  1.       ]))
    np.testing.assert_array_almost_equal(space.nodes[1],
                                         np.array([-1.00000000e+00,  1.23259516e-32,  1.00000000e+00]))

    quad_nodes, quad_weights, quad_weights_ravel = space.___PRIVATE_do_evaluate_quadrature___((3, 2))

    np.testing.assert_array_almost_equal(
        quad_nodes[0], np.array([-0.86113631, -0.33998104,  0.33998104,  0.86113631]))
    np.testing.assert_array_almost_equal(
        quad_nodes[1], np.array([-0.77459667,  0.        ,  0.77459667]))

    np.testing.assert_array_almost_equal(
        quad_weights[0], np.array([0.34785485, 0.65214515, 0.65214515, 0.34785485]))
    np.testing.assert_array_almost_equal(
        quad_weights[1], np.array([0.55555556, 0.88888889, 0.55555556]))

    np.testing.assert_array_almost_equal(
        quad_weights_ravel, np.array([0.19325269, 0.36230286, 0.36230286, 0.19325269, 0.30920431,
                                      0.57968458, 0.57968458, 0.30920431, 0.19325269, 0.36230286,
                                      0.36230286, 0.19325269]))

    assert space.Kronecker == (True, [True, True])
    assert space.IS_Kronecker
    assert space.standard_properties.stamp == '2dCSCG|structured|space'

    return 1






if __name__ == '__main__':
    # mpiexec python _2dCSCG\TESTS\unittest_space.py
    test_Space_NO1_polynomial_space()