# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

To run all test, do
::

    $ mpiexec python _3dCSCG\TESTS\mpi_unittest.py

This calls all cores. To call particular number of cores, say ``3`` cores, do
::

    $ mpiexec -n 3 python _3dCSCG\TESTS\mpi_unittest.py

It is always suggested testing the library multiple time with different numbers of cores.

"""
import sys
if './' not in sys.path: sys.path.append('./')
passed_3dCSCG_tests = 0

from _3dCSCG.TESTS.unittest_forms import *
from _3dCSCG.TESTS.unittest_fields import *
from _3dCSCG.TESTS.unittest_spaces import *
from _3dCSCG.TESTS.unittest_mesh import *
from _3dCSCG.TESTS.unittest_Naive_numbering import *
from _3dCSCG.TESTS.unittest_APP import *
from _3dCSCG.TESTS.unittest_ADF import *
from _3dCSCG.TESTS.unittest_trace import *
from _3dCSCG.TESTS.unittest_edge_forms import *

t_3dCSCG_start = MPI.Wtime()

if rAnk == mAster_rank: print(f"\n [_3dCSCG] tests start...\n")


passed_3dCSCG_tests += test_ADF_NO1_general_tests_standard_forms()
passed_3dCSCG_tests += test_ADF_NO2_general_tests_trace_forms()
passed_3dCSCG_tests += test_ADF_NO3_coboundary()

passed_3dCSCG_tests += test_APP_NO1_icpsNS_no_hybrid_manu_conserving()
passed_3dCSCG_tests += test_APP_NO2_icpsNS_no_hybrid_TGV()

passed_3dCSCG_tests += test_Form_NO0_3dCSCG_Field_numerical()
passed_3dCSCG_tests += test_Form_NO1_3dCSCG_VectorField()
passed_3dCSCG_tests += test_Form_NO2_3dCSCG_ScalarField()
passed_3dCSCG_tests += test_Form_NO3_3dCSCG_TensorField()

passed_3dCSCG_tests += test_Form_NO1_discretization_and_reconstruction()
passed_3dCSCG_tests += test_Form_NO1a_discretization_and_reconstruction()
passed_3dCSCG_tests += test_Form_NO1b_trace_form_Rd_and_Rc()
passed_3dCSCG_tests += test_Form_NO2_mass_matrix()
passed_3dCSCG_tests += test_Form_NO3_incidence_matrices()
passed_3dCSCG_tests += test_Form_NO4_cross_product_1()
passed_3dCSCG_tests += test_Form_NO5_cross_product_2()
passed_3dCSCG_tests += test_Form_NO6_resemble()
passed_3dCSCG_tests += test_Form_No7_with_other_element_numbering_AUTO()
passed_3dCSCG_tests += test_Form_No8_edge_forms()
passed_3dCSCG_tests += test_Form_No9_node_forms()


passed_3dCSCG_tests += test_Space_NO1_basis_functions_mapping_test()

passed_3dCSCG_tests += test_Mesh_NO0_element_division_and_numbering_quality()
passed_3dCSCG_tests += test_Mesh_NO1_mesh_general()
passed_3dCSCG_tests += test_Mesh_NO2_trace_elements()
passed_3dCSCG_tests += test_Mesh_NO2a_trace_elements_CT()
passed_3dCSCG_tests += test_Mesh_NO3_elements_CT()
passed_3dCSCG_tests += test_Mesh_NO4_elements_CT_QUAD()
passed_3dCSCG_tests += test_Mesh_NO5_mesh_trace_topology()
passed_3dCSCG_tests += test_Mesh_NO5a_mesh_trace_CT()
passed_3dCSCG_tests += test_Mesh_NO6_transfinite()
passed_3dCSCG_tests += test_Mesh_NO7_boundaries()
passed_3dCSCG_tests += test_Mesh_NO8_Mesh_SubGeometry_perpendicular_slice_object()
passed_3dCSCG_tests += test_Mesh_NO9_edge_node_mesh()
passed_3dCSCG_tests += test_Form_No10_standard_form_dofs()

passed_3dCSCG_tests += test_Naive_Numbering_NO1_0form()
passed_3dCSCG_tests += test_Naive_Numbering_NO2_1form()
passed_3dCSCG_tests += test_Naive_Numbering_NO3_2form()
passed_3dCSCG_tests += test_Naive_Numbering_NO4_2trace()
passed_3dCSCG_tests += test_Naive_Numbering_NO6_1trace()
passed_3dCSCG_tests += test_Naive_Numbering_NO5_0trace()


passed_3dCSCG_tests += test_trace_NO__general_tests()
passed_3dCSCG_tests += test_trace_NO0_trace_0_form_Rd_and_Rc()
passed_3dCSCG_tests += test_trace_NO1_trace_1_form_Rd_and_Rc()
passed_3dCSCG_tests += test_trace_NO2_trace_2_form_Rd_and_Rc()
passed_3dCSCG_tests += test_trace_NO3_trace_matrices()

passed_3dCSCG_tests += test_edge_forms_No0_save_read()
passed_3dCSCG_tests += test_edge_forms_No1_0edge_Rd_and_Rc()


if rAnk == mAster_rank:
    print("\n<{}> _3dCSCG tests passed; cost {:.3f} seconds.\n".format(
        passed_3dCSCG_tests, MPI.Wtime()-t_3dCSCG_start))