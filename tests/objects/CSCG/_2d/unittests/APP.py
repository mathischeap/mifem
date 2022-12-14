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
from tests.objects.CSCG._2d.unittests.auxiliaries.Poisson_essential_BC import PoissonSolver1
from tests.objects.CSCG._2d.unittests.auxiliaries.scalar_Laplace_essential_BC_iterative_solver import \
    scalar_Laplace_solver_iterative_solver

from tests.objects.CSCG._2d.unittests.auxiliaries.Euler_shear_layer_rollup_direct import \
    Euler_shear_layer_rollup_direct_test

from components.miscellaneous.mios import isfile, remove



def test_APP_NO2_scalar_Laplace_essential_BC_iterative_solver():
    """"""
    if RANK == MASTER_RANK:
        print(">>> [test_APP_NO2_scalar_Laplace_essential_BC_iterative_solver] ...... ", flush=True)

    u_error_L2, p_error_L2 = scalar_Laplace_solver_iterative_solver(0.0, 4, 3, 4, 5)

    assert u_error_L2 < 0.005
    assert p_error_L2 < 0.0007

    return 1



def test_APP_NO3_Euler_ShearLayerRollup_Direct_test():
    """"""
    if RANK == MASTER_RANK:
        print(">>> [test_APP_NO3_Euler_ShearLayerRollup_Direct_test] ...... ", flush=True)
    K = 4 # K * K elements (uniform)
    N = 2  # polynomial degree
    dt = 0.2
    t = 1
    image_folder = './APP_test_No3_images_direct'
    RDF_filename = 'shear_layer_rollup_direct_test'

    if isfile(RDF_filename + '.csv'): remove(RDF_filename + '.csv')

    SI = Euler_shear_layer_rollup_direct_test(K, N, dt, t, image_folder, RDF_filename)

    if RANK == MASTER_RANK:
        os.remove(image_folder + '/video.avi')
        os.rmdir(image_folder)
        os.remove(RDF_filename + '.csv')

        data = SI.RDF.to_numpy()

        np.testing.assert_array_almost_equal(data[-1,:],
                np.array([1.000000e+00, 2.000000e-01, 4.239597e+01, 1.684394e+01,
                          1.776357e-15, 3.052409e-14]), decimal=5)

    return 1



def test_APP_NO4_Poisson_hMSEM_test_1():
    """"""
    if RANK == MASTER_RANK:
        print(">>> [test_APP_NO4_Poisson_hMSEM_test_1] ...... ", flush=True)

    p_error_L2, u_error_L2 = PoissonSolver1(0.125, 8, 7, 3, 4)

    assert p_error_L2 < 0.0008
    assert u_error_L2 < 0.002

    return 1



if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/__tests__/unittests/APP.py
    test_APP_NO3_Euler_ShearLayerRollup_Direct_test()
    # test_APP_NO2_scalar_Laplace_essential_BC_iterative_solver()

