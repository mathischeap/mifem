# -*- coding: utf-8 -*-
"""
mpiexec python objects/CSCG/_2d/__tests__/unittests/main.py
mpiexec -n 4 python objects/CSCG/_2d/__tests__/unittests/main.py

"""

import sys
if './' not in sys.path: sys.path.append('./')
passed_2dCSCG_tests = 0

from objects.CSCG._2d.__tests__.unittests.standard_forms.general import *
from objects.CSCG._2d.__tests__.unittests.standard_forms.dofs import test_standard_forms_DOFS
from objects.CSCG._2d.__tests__.unittests.mesh import *
from objects.CSCG._2d.__tests__.unittests.spaces import *
from objects.CSCG._2d.__tests__.unittests.APP import *
from objects.CSCG._2d.__tests__.unittests.fields import *

t_3dCSCG_start = MPI.Wtime()

if rAnk == mAster_rank:
    print(f"\n [_2dCSCG] tests start...\n")


passed_2dCSCG_tests += test_Form_NO1_coboundary()
passed_2dCSCG_tests += test_Form_NO2_Naive_numbering()
passed_2dCSCG_tests += test_Form_NO3_mass_matrices()
passed_2dCSCG_tests += test_Form_NO4_Hodge()
passed_2dCSCG_tests += test_Form_NO5_cross_product()
passed_2dCSCG_tests += test_Form_NO6_reconstruction_matrices()
passed_2dCSCG_tests += test_Form_NO7_weak_curl()
passed_2dCSCG_tests += test_standard_forms_DOFS()

passed_2dCSCG_tests += test_Mesh_NO1_mesh_topology()
passed_2dCSCG_tests += test_Mesh_NO2_mesh_coordinate_transformation()
passed_2dCSCG_tests += test_Mesh_NO3_mesh_coordinate_transformation_QUAD()
passed_2dCSCG_tests += test_Mesh_NO4_mesh_trace_topology()

passed_2dCSCG_tests += test_Space_NO1_polynomial_space()

passed_2dCSCG_tests += test_APP_NO1_scalar_Laplace_essential_BC()
passed_2dCSCG_tests += test_APP_NO2_scalar_Laplace_essential_BC_iterative_solver()
passed_2dCSCG_tests += test_APP_NO3_Euler_ShearLayerRollup_Direct_test()

passed_2dCSCG_tests += test_Fields_NO1_vector()
passed_2dCSCG_tests += test_Fields_NO2_scalar()




if rAnk == mAster_rank:
    print("\n<{}> _2dCSCG tests passed; cost {:.3f} seconds.\n".format(
        passed_2dCSCG_tests, MPI.Wtime()-t_3dCSCG_start))