# -*- coding: utf-8 -*-
"""
Here we put all unittests in mpi environment here.

mpiexec -n 4 python __tests__/test_all.py

mpiexec -n 4 python __tests__/unittests/main.py

mpiexec -n 4 python objects/CSCG/_2d/__tests__/unittests/main.py

mpiexec -n 4 python objects/CSCG/_3d/__tests__/unittests/main.py

"""

import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *

if rAnk == mAster_rank:
    from screws.miscellaneous.timer import count_files_and_lines
    files, lines = count_files_and_lines('./')

passed_2dCSCG_tests = 0
passed_3dCSCG_tests = 0
passed_GLOBAL_tests = 0
t_global_start = MPI.Wtime()




from objects.CSCG._2d.__tests__.unittests.main import passed_2dCSCG_tests # comment to skip these tests.


from objects.CSCG._3d.__tests__.unittests.main import passed_3dCSCG_tests # comment to skip these tests.


from __tests__.unittests.main import passed_GLOBAL_tests # comment to skip these tests.






total_Tests = passed_2dCSCG_tests + passed_3dCSCG_tests + passed_GLOBAL_tests

if rAnk == mAster_rank:
    print("\n<{}> total tests passed; cost {:.3f} seconds.\n".format(
        total_Tests, MPI.Wtime()-t_global_start))