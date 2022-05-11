# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, the Netherlands

"""
import os

import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from objects.CSCG._2d.__tests__.unittests.auxiliaries.scalar_Laplace_essential_BC import scalar_Laplace_solver
from objects.CSCG._2d.__tests__.unittests.auxiliaries.scalar_Laplace_essential_BC_iterative_solver import \
    scalar_Laplace_solver_iterative_solver

from objects.CSCG._2d.__tests__.unittests.auxiliaries.Euler_shear_layer_rollup_direct import \
    Euler_shear_layer_rollup_direct_test



def test_APP_NO1_scalar_Laplace_essential_BC():
    """"""
    if rAnk == mAster_rank:
        print(">>> [test_APP_NO1_scalar_Laplace_essential_BC] ...... ", flush=True)

    u_error_L2, p_error_L2 = scalar_Laplace_solver(0.15, 4, 3, 4, 5)

    assert u_error_L2 < 3e-2
    assert p_error_L2 < 8e-3

    return 1


def test_APP_NO2_scalar_Laplace_essential_BC_iterative_solver():
    """"""
    if rAnk == mAster_rank:
        print(">>> [test_APP_NO2_scalar_Laplace_essential_BC_iterative_solver] ...... ", flush=True)

    u_error_L2, p_error_L2 = scalar_Laplace_solver_iterative_solver(0.0, 4, 3, 4, 5)

    assert u_error_L2 < 0.005
    assert p_error_L2 < 0.0007

    return 1



def test_APP_NO3_Euler_ShearLayerRollup_Direct_test():
    """"""
    if rAnk == mAster_rank:
        print(">>> [test_APP_NO3_Euler_ShearLayerRollup_Direct_test] ...... ", flush=True)
    K = 4 # K * K elements (uniform)
    N = 2  # polynomial degree
    dt = 0.2
    t = 1
    image_folder = './APP_test_No3_images_direct'
    RDF_filename = 'shear_layer_rollup_direct_test'

    SI = Euler_shear_layer_rollup_direct_test(K, N, dt, t, image_folder, RDF_filename)

    if rAnk == mAster_rank:
        os.remove(image_folder + '/video.avi')
        os.rmdir(image_folder)
        os.remove(RDF_filename + '.csv')

        data = SI.RDF.to_numpy()

        np.testing.assert_array_almost_equal(data[-1,:],
                np.array([1.00000000e+00, 2.00000000e-01, 4.23959684e+01, 1.68439419e+01,
                          0.00000000e+00, 1.65772729e-14]),)

    return 1








if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/__tests__/unittests/APP.py
    test_APP_NO3_Euler_ShearLayerRollup_Direct_test()
    # test_APP_NO2_scalar_Laplace_essential_BC_iterative_solver()

