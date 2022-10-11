# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

To run all test with given number of threads, do

>> mpiexec -n 6 python __tests__/unittests/main.py

It is always suggested testing the library multiple time with different numbers of cores.

"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *

passed_GLOBAL_tests = 0

if rAnk == mAster_rank:
    print(f"\n [Global] tests start...\n")

from __tests__.unittests.mifem import *
from __tests__.unittests.tools_ import *
from __tests__.unittests.screws_ import *
from __tests__.unittests.linear_solvers import *
from __tests__.unittests.EWC.column_vector import test_LinearAlgebra_EWC_No0_ColumnVector
from __tests__.unittests.EWC.operators import test_LinearAlgebra_EWC_No1_Operators
from __tests__.unittests.EWC.adjust import test_EWC_adjust

from __tests__.unittests.gathering_matrix.find import test_GatheringMatrix_find
from __tests__.unittests.gathering_matrix.chain_gathering_matrix import Test_ChainGM_sequent_chain_method

from __tests__.unittests.linear_system.strong_BC_of_Poisson import \
    test_applying_strong_BC_for_Poisson_problem
from __tests__.unittests.linear_system.strong_BC_of_Poisson_non_trivial_BC import \
    test_applying_strong_BC_for_Poisson_problem_NT
from __tests__.unittests.MultiDimMatrix.cross_product_MDM_test import test_MDM_sf_CrossProduct
from __tests__.unittests.nonlinear_solver.regular_Newton_Raphson import test_Regular_Newton_Raphson
from __tests__.unittests.Parallel_Matrix3dInput_Runner.WellTest import WellTest_ParallelMatrix3dInputRunner

from __tests__.unittests.VTK.unstructuredGridToVTK import TEST_save_CSCG_objects_to_unstructured_VTK_file
from __tests__.unittests.VTK.gridToVTK import TEST_save_CSCG_objects_to_structured_VTK_file

from __tests__.unittests.gathering_matrix.customize_sequent import CustomizeSequent



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
passed_GLOBAL_tests += test_LinearSolver_No3_direct()

passed_GLOBAL_tests += test_LinearAlgebra_EWC_No0_ColumnVector()
passed_GLOBAL_tests += test_LinearAlgebra_EWC_No1_Operators()
passed_GLOBAL_tests += test_EWC_adjust()

passed_GLOBAL_tests += test_GatheringMatrix_find()
passed_GLOBAL_tests += Test_ChainGM_sequent_chain_method()()
passed_GLOBAL_tests += CustomizeSequent()()
passed_GLOBAL_tests += test_applying_strong_BC_for_Poisson_problem()
passed_GLOBAL_tests += test_applying_strong_BC_for_Poisson_problem_NT()

passed_GLOBAL_tests += test_MDM_sf_CrossProduct()
passed_GLOBAL_tests += test_Regular_Newton_Raphson()
passed_GLOBAL_tests += WellTest_ParallelMatrix3dInputRunner()

passed_GLOBAL_tests += TEST_save_CSCG_objects_to_structured_VTK_file()
passed_GLOBAL_tests += TEST_save_CSCG_objects_to_unstructured_VTK_file()



cOmm.barrier()
if rAnk == mAster_rank:
    print("\n<{}> Global tests passed; cost {:.3f} seconds.\n".format(
        passed_GLOBAL_tests, MPI.Wtime()-t_global_start))