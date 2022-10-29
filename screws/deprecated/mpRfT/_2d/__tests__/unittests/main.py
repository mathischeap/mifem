# -*- coding: utf-8 -*-
"""

mpiexec -n 4 python objects/mpRfT/_2d/__tests__/unittests/main.py

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 3:31 PM
"""
import sys


if './' not in sys.path: sys.path.append('./')
from root.config.main import MPI, RANK, MASTER_RANK


passed_mpRfT2_tests = 0
start_time = MPI.Wtime()
if RANK == MASTER_RANK: print(f"\n [nCSCG_RF2] tests start...\n")


from objects.mpRfT._2d.__tests__.unittests.save_read import test_mpRfT2_save_read
from objects.mpRfT._2d.__tests__.unittests.numbering import test_mpRfT2_ir_numbering


passed_mpRfT2_tests += test_mpRfT2_save_read()
passed_mpRfT2_tests += test_mpRfT2_ir_numbering()


if RANK == MASTER_RANK:
    print("\n<{}> mpRfT2 tests passed; cost {:.3f} seconds.\n".format(
        passed_mpRfT2_tests, MPI.Wtime()-start_time))