# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config import *
from _2dCSCG.TESTS.__unittest_scripts__.scalar_Laplace_essential_BC import scalar_Laplace_solver


def test_APP_NO1_scalar_Laplace_essential_BC():
    """"""
    if rAnk == mAster_rank:
        print(">>> [test_APP_NO1_scalar_Laplace_essential_BC] ...... ", flush=True)

    u_error_L2, p_error_L2 = scalar_Laplace_solver(0.15, 4, 3, 4, 5)

    assert u_error_L2 < 3e-2
    assert p_error_L2 < 8e-3

    return 1


if __name__ == '__main__':
    # mpiexec python _2dCSCG\TESTS\unittest_app.py
    test_APP_NO1_scalar_Laplace_essential_BC()

