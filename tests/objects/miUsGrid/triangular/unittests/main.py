# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

To run all test, do

mpiexec python tests/objects/miUsGrid/triangular/unittests/main.py

This calls all cores. To call particular number of cores, say ``3`` cores, do

mpiexec -n 3 python tests/objects/miUsGrid/triangular/unittests/main.py

It is always suggested testing the library multiple time with different numbers of cores.

"""
import sys
if './' not in sys.path: sys.path.append('./')
passed_miUSGridTriangle_tests = 0

from root.config.main import MPI, RANK, MASTER_RANK, COMM

t_start = MPI.Wtime()

from tests.objects.miUsGrid.triangular.unittests.standardForms.incidence_matrices import \
    miUsGrid_Triangle_Incidence_matrices
from tests.objects.miUsGrid.triangular.unittests.standardForms.convergence_test.test import \
    miUsGrid_TriangleMesh_ConvergenceTest
from tests.objects.miUsGrid.triangular.unittests.standardForms.dofs_topology import \
    Test_dofs_topology_S1F
from tests.objects.miUsGrid.triangular.unittests.standardForms.reconstruction_matrices import \
    miUsGrid_Triangles_ReconstructionMatrices
from tests.objects.miUsGrid.triangular.unittests.standardForms.Poisson.test import \
    miUsGrid_Triangle_Poisson
from tests.objects.miUsGrid.triangular.unittests.standardForms.mass_matrices import \
    miUsGrid_Triangles_MassMatrices
from tests.objects.miUsGrid.triangular.unittests.MSEM_Stokes_test import miUsTriangleTest_MSEM_STOKES

passed_miUSGridTriangle_tests += miUsGrid_Triangle_Incidence_matrices()()
passed_miUSGridTriangle_tests += miUsGrid_TriangleMesh_ConvergenceTest()()
passed_miUSGridTriangle_tests += Test_dofs_topology_S1F()()
passed_miUSGridTriangle_tests += miUsGrid_Triangles_ReconstructionMatrices()()
passed_miUSGridTriangle_tests += miUsGrid_Triangle_Poisson()()
passed_miUSGridTriangle_tests += miUsGrid_Triangles_MassMatrices()()
passed_miUSGridTriangle_tests += miUsTriangleTest_MSEM_STOKES()

COMM.barrier()
if RANK == MASTER_RANK:
    print("\n<{}> miUSGrid Triangle tests passed; cost {:.3f} seconds.\n".format(
        passed_miUSGridTriangle_tests, MPI.Wtime()-t_start))