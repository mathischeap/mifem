# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, the Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from _2dCSCG.tests.unittests.auxiliaries.scalar_Laplace_essential_BC import scalar_Laplace_solver
from _2dCSCG.tests.unittests.auxiliaries.scalar_Laplace_essential_BC_iterative_solver import \
    scalar_Laplace_solver_iterative_solver





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



if __name__ == '__main__':
    # mpiexec -n 4 python _2dCSCG\tests\unittests\APP.py
    test_APP_NO1_scalar_Laplace_essential_BC()
    test_APP_NO2_scalar_Laplace_essential_BC_iterative_solver()

