# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

To run all test with given number of threads, do

>> mpiexec -n 6 python tests/unittests/main.py

It is always suggested testing the library multiple time with different numbers of cores.

"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *

if rAnk == mAster_rank:
    from screws.miscellaneous.timer import count_files_and_lines
    files, lines = count_files_and_lines('./')

t_start = MPI.Wtime()

passed_2dCSCG_tests = 0
passed_3dCSCG_tests = 0
passed_GLOBAL_tests = 0

from _2dCSCG.tests.unittests.main import * # comment this line to skip these tests.
from _3dCSCG.tests.unittests.main import * # comment this line to skip these tests.

if rAnk == mAster_rank:
    print(f"\n [Global] tests start...\n")



from tests.unittests.mifem import *
from tests.unittests.tools_ import *
from tests.unittests.screws_ import *
from tests.unittests.linear_solvers import *
from tests.unittests.linear_algebra import *

t_global_start = MPI.Wtime()

passed_GLOBAL_tests += test_mifem_NO1_2dCSCG_save_read()
passed_GLOBAL_tests += test_mifem_NO2_3dCSCG_save_read()

passed_GLOBAL_tests += test_TOOLS_NO1_iterator()

passed_GLOBAL_tests += test_TOOLS_NO2_0_linear_algebra_gmres0_solver()
passed_GLOBAL_tests += test_TOOLS_NO2_1_linear_algebra_gmres1_solver()
passed_GLOBAL_tests += test_TOOLS_NO2_2_linear_algebra_gmres2_solver()
passed_GLOBAL_tests += test_TOOLS_NO2_3_linear_algebra_serial_scipy_sparse_solver()
passed_GLOBAL_tests += test_TOOLS_NO2_4_linear_algebra_serial_spsolve_solver()

passed_GLOBAL_tests += test_TOOLS_NO3_GlobalMatrix_GlobalVector_operators_test()
passed_GLOBAL_tests += test_TOOLS_NO4_GlobalMatrix_GlobalVector_operators_test()
passed_GLOBAL_tests += test_TOOLS_NO5_DistributedVector_operators_test()
passed_GLOBAL_tests += test_TOOLS_NO6_send_GM_in_parts_test()
passed_GLOBAL_tests += test_TOOLS_NO7_linear_algebra_EWC_test()
passed_GLOBAL_tests += test_TOOLS_NO8_GlobalMatrix_dot_product_test()
passed_GLOBAL_tests += test_TOOLS_NO9_test_Chained_Gathering_Matrix()
passed_GLOBAL_tests += test_TOOLS_NO10_test_EWC_SparseMatrix_Customize()
passed_GLOBAL_tests += test_TOOLS_NO11_test_ParallelMatrix3dInputRunner()
passed_GLOBAL_tests += test_TOOLS_NO12_EWC_assembling_test()
passed_GLOBAL_tests += test_TOOLS_NO13_EWC_Customize_CSCG_partial_dofs()
passed_GLOBAL_tests += test_TOOLS_NO14_partial_cochain_with_3dCSCG_form_BC()
passed_GLOBAL_tests += test_TOOLS_NO15_linear_system_apply_BC()

passed_GLOBAL_tests += test_SCREWS_NO1_3d_functions()
passed_GLOBAL_tests += test_SCREWS_NO2_sending_an_email_to_admin()

passed_GLOBAL_tests += test_LinearSolver_No0_GMRES()
passed_GLOBAL_tests += test_LinearSolver_No1_BiCGSTAB()
passed_GLOBAL_tests += test_LinearSolver_No2_LooseGMRES()

passed_GLOBAL_tests += test_LinearAlgebra_No0_EWC_ColumnVector()


total_Tests = passed_2dCSCG_tests + passed_3dCSCG_tests + passed_GLOBAL_tests
cOmm.barrier()
if rAnk == mAster_rank:
    print("\n<{}> Global tests passed; cost {:.3f} seconds.\n".format(
        passed_GLOBAL_tests, MPI.Wtime()-t_global_start))
    print("\nIN TOTAL, <{}> tests passed; cost {:.3f} seconds.\n".format(
        total_Tests, MPI.Wtime()-t_start))