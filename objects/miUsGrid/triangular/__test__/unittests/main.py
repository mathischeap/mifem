# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

To run all test, do

mpiexec python objects/miUsGrid/triangular/__test__/unittests/main.py

This calls all cores. To call particular number of cores, say ``3`` cores, do

mpiexec -n 3 python objects/miUsGrid/triangular/__test__/unittests/main.py

It is always suggested testing the library multiple time with different numbers of cores.

"""
import sys
if './' not in sys.path: sys.path.append('./')
passed_miUSGridTriangle_tests = 0

from root.config.main import MPI, rAnk, mAster_rank, cOmm

t_start = MPI.Wtime()


from objects.miUsGrid.triangular.__test__.unittests.standard_forms.incidence_matrices import \
    miUsGrid_Triangle_Incidence_matrices
from objects.miUsGrid.triangular.__test__.unittests.standard_forms.convergence_test.test import \
    miUsGrid_TriangleMesh_ConvergenceTest
from objects.miUsGrid.triangular.__test__.unittests.standard_forms.dofs_topology import \
    Test_dofs_topology_S1F
from objects.miUsGrid.triangular.__test__.unittests.standard_forms.reconstruction_matrices import \
    miUsGrid_Triangles_ReconstructionMatrices
from objects.miUsGrid.triangular.__test__.unittests.standard_forms.Poisson.test import \
    miUsGrid_Triangle_Poisson
from objects.miUsGrid.triangular.__test__.unittests.standard_forms.mass_matrices import \
    miUsGrid_Triangles_MassMatrices
from objects.miUsGrid.triangular.__test__.unittests.MSEM_Stokes_test import miUsTriangleTest_MSEM_STOKES


passed_miUSGridTriangle_tests += miUsGrid_Triangle_Incidence_matrices()()
passed_miUSGridTriangle_tests += miUsGrid_TriangleMesh_ConvergenceTest()()
passed_miUSGridTriangle_tests += Test_dofs_topology_S1F()()
passed_miUSGridTriangle_tests += miUsGrid_Triangles_ReconstructionMatrices()()
passed_miUSGridTriangle_tests += miUsGrid_Triangle_Poisson()()
passed_miUSGridTriangle_tests += miUsGrid_Triangles_MassMatrices()()
passed_miUSGridTriangle_tests += miUsTriangleTest_MSEM_STOKES()


cOmm.barrier()
if rAnk == mAster_rank:
    print("\n<{}> miUSGrid Triangle tests passed; cost {:.3f} seconds.\n".format(
        passed_miUSGridTriangle_tests, MPI.Wtime()-t_start))