# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

To run all 3d test, do

>> mpiexec python objects/CSCG/_3d/__tests__/unittests/main.py

This calls all cores. To call particular number of cores, say ``3`` cores, do

>> mpiexec -n 3 python objects/CSCG/_3d/__tests__/unittests/main.py

It is always suggested testing the library multiple time with different numbers of cores.

"""

import sys
if './' not in sys.path: sys.path.append('./')
passed_3dCSCG_tests = 0

from tests.objects.CSCG._3d.unittests.standard_forms.general import *
from tests.objects.CSCG._3d.unittests.standard_forms.reconstruct_2_DF import *
from tests.objects.CSCG._3d.unittests.fields import *
from tests.objects.CSCG._3d.unittests.spaces import *
from tests.objects.CSCG._3d.unittests.mesh import *
from tests.objects.CSCG._3d.unittests.Naive_numbering import *
from tests.objects.CSCG._3d.unittests.APP import *
from tests.objects.CSCG._3d.unittests.ADF import *
from tests.objects.CSCG._3d.unittests.trace_forms import *
from tests.objects.CSCG._3d.unittests.edge_forms import *
from tests.objects.CSCG._3d.unittests.hybrid.edge0 import test_hybridization_of_standard_0_form
from tests.objects.CSCG._3d.unittests.hybrid.edge1 import test_hybridization_of_standard_1_form
from tests.objects.CSCG._3d.unittests.hybrid.trace2 import test_hybridization_trace2
from tests.objects.CSCG._3d.unittests.trace_selective_matrices import test_trace_and_selective_matrices
from tests.objects.CSCG._3d.unittests.Stokes_flow.MSEM_trivial_BC import test_Stokes_MSEM_trivial_BC
from tests.objects.CSCG._3d.unittests.Stokes_flow.hdMSEM_direct_solver import test_Stokes_hdMSEM
from tests.objects.CSCG._3d.unittests.Stokes_flow.hdMSEM_Schur_direct_3 import test_Stokes_hdMSEM_Schur_Rank2Solver

from tests.objects.CSCG._3d.unittests.Poisson.hdMSEM_inner import test_hdMSEM_Poisson_Inner

from tests.objects.CSCG._3d.unittests.TISE.hdMSEM_inner import test_hdMSEM_Schrodinger_Inner
from tests.objects.CSCG._3d.unittests.TISE.hdMSEM_inner_cuboid import test_hdMSEM_Schrodinger_Inner_Cuboid

t_3dCSCG_start = MPI.Wtime()

if RANK == MASTER_RANK: print(f"\n [_3dCSCG] tests start...\n")

passed_3dCSCG_tests += test_ADF_NO1_general_tests_standard_forms()
passed_3dCSCG_tests += test_ADF_NO2_general_tests_trace_forms()
passed_3dCSCG_tests += test_ADF_NO3_coboundary()

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
passed_3dCSCG_tests += test_Form_NOx_cross_product_3()
passed_3dCSCG_tests += test_Form_NOx1_cross_product_4()
passed_3dCSCG_tests += test_Form_NOx2_F_dot_G_times_H()
passed_3dCSCG_tests += test_Form_No7_with_other_element_numbering_AUTO()
passed_3dCSCG_tests += test_Form_No10_standard_form_dofs()
passed_3dCSCG_tests += test_Form_No11_reconstruction_matrices()
passed_3dCSCG_tests += test_Form_NO12_weak_curl()

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

passed_3dCSCG_tests += test_hybridization_of_standard_0_form()
passed_3dCSCG_tests += test_hybridization_of_standard_1_form()
passed_3dCSCG_tests += test_hybridization_trace2()
passed_3dCSCG_tests += test_trace_and_selective_matrices()
passed_3dCSCG_tests += test_reconstruct_DF()

from components.miscellaneous.mirand import randint
a = randint(0, 3) # below tests are slow, we do not run them at every single time.
if a == 1:
    passed_3dCSCG_tests += test_Stokes_MSEM_trivial_BC()
    passed_3dCSCG_tests += test_Stokes_hdMSEM()
    passed_3dCSCG_tests += test_Stokes_hdMSEM_Schur_Rank2Solver()
elif a == 2:
    passed_3dCSCG_tests += test_hdMSEM_Poisson_Inner()
elif a == 3:
    passed_3dCSCG_tests += test_hdMSEM_Schrodinger_Inner()
    passed_3dCSCG_tests += test_hdMSEM_Schrodinger_Inner_Cuboid()
else:
    pass


if RANK == MASTER_RANK:
    print("\n<{}> _3dCSCG tests passed; cost {:.3f} seconds.\n".format(
        passed_3dCSCG_tests, MPI.Wtime()-t_3dCSCG_start))