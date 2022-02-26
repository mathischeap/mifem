"""
# -*- coding: utf-8 -*-

>> mpiexec python _2dCSCG/tests/unittests/main.py
>> mpiexec -n 3 python _2dCSCG/tests/unittests/main.py

"""

import sys
if './' not in sys.path: sys.path.append('./')
passed_2dCSCG_tests = 0

from _2dCSCG.tests.unittests.standard_forms import *
from _2dCSCG.tests.unittests.mesh import *
from _2dCSCG.tests.unittests.spaces import *
from _2dCSCG.tests.unittests.APP import *

t_3dCSCG_start = MPI.Wtime()

if rAnk == mAster_rank:
    print(f"\n [_2dCSCG] tests start...\n")


passed_2dCSCG_tests += test_Form_NO1_coboundary()
passed_2dCSCG_tests += test_Form_NO2_Naive_numbering()
passed_2dCSCG_tests += test_Form_NO3_mass_matrices()
passed_2dCSCG_tests += test_Form_NO4_Hodge()


passed_2dCSCG_tests += test_Mesh_NO1_mesh_topology()
passed_2dCSCG_tests += test_Mesh_NO2_mesh_coordinate_transformation()
passed_2dCSCG_tests += test_Mesh_NO3_mesh_coordinate_transformation_QUAD()
passed_2dCSCG_tests += test_Mesh_NO4_mesh_trace_topology()


passed_2dCSCG_tests += test_Space_NO1_polynomial_space()


passed_2dCSCG_tests += test_APP_NO1_scalar_Laplace_essential_BC()


if rAnk == mAster_rank:
    print("\n<{}> _2dCSCG tests passed; cost {:.3f} seconds.\n".format(
        passed_2dCSCG_tests, MPI.Wtime()-t_3dCSCG_start))