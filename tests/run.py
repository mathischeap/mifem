# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

mpiexec -n 4 python tests/run.py

mpiexec -n 4 python tests/tools/main.py

mpiexec -n 4 python tests/objects/CSCG/_2d/unittests/main.py

mpiexec -n 4 python tests/objects/CSCG/_3d/unittests/main.py

mpiexec -n 4 python tests/objects/miUsGrid/triangular/unittests/main.py

"""

import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import RANK, MASTER_RANK, MPI

if RANK == MASTER_RANK:
    from screws.miscellaneous.timer import count_files_and_lines
    files, lines = count_files_and_lines('./')


passed_2dCSCG_tests = 0           # do not comment this
passed_3dCSCG_tests = 0           # do not comment this
passed_GLOBAL_tests = 0           # do not comment this
passed_miUSGridTriangle_tests = 0 # do not comment this

t_global_start = MPI.Wtime()



from tests.tools.main import passed_GLOBAL_tests # comment to skip these tests.


from tests.objects.CSCG._2d.unittests.main import passed_2dCSCG_tests # comment to skip these tests.


from tests.objects.CSCG._3d.unittests.main import passed_3dCSCG_tests # comment to skip these tests.


from tests.objects.miUsGrid.triangular.unittests.main import passed_miUSGridTriangle_tests


total_Tests = passed_2dCSCG_tests + \
              passed_3dCSCG_tests + \
              passed_GLOBAL_tests + \
              passed_miUSGridTriangle_tests



if RANK == MASTER_RANK:
    print("\n<{}> total tests passed; cost {:.3f} seconds.\n".format(
        total_Tests, MPI.Wtime()-t_global_start))