# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

To run all test, do

>> mpiexec python objects/miUsGrid/triangular/__test__/unittests/main.py

This calls all cores. To call particular number of cores, say ``3`` cores, do

>> mpiexec -n 3 python objects/miUsGrid/triangular/__test__/unittests/main.py

It is always suggested testing the library multiple time with different numbers of cores.

"""
import sys
if './' not in sys.path: sys.path.append('./')
passed_miUSGridTriangle_tests = 0

from root.config.main import MPI, rAnk, mAster_rank

t_start = MPI.Wtime()


from objects.miUsGrid.triangular.__test__.unittests.standard_forms.incidence_matrices import \
    miUsGrid_Triangle_Incidence_matrices
from objects.miUsGrid.triangular.__test__.unittests.standard_forms.convergence_test.test import \
    miUsGrid_TriangleMesh_ConvergenceTest
from objects.miUsGrid.triangular.__test__.unittests.standard_forms.numbering import \
    miUsGrid_Triangle_StandardFormNumbering

passed_miUSGridTriangle_tests += miUsGrid_Triangle_Incidence_matrices()()
passed_miUSGridTriangle_tests += miUsGrid_TriangleMesh_ConvergenceTest()()
passed_miUSGridTriangle_tests += miUsGrid_Triangle_StandardFormNumbering()()

if rAnk == mAster_rank:
    print("\n<{}> miUSGrid Triangle tests passed; cost {:.3f} seconds.\n".format(
        passed_miUSGridTriangle_tests, MPI.Wtime()-t_start))