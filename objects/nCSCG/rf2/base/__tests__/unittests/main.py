# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 3:14 PM

mpiexec -n 4 python objects/nCSCG/rf2/base/__tests__/unittests/main.py

"""
import sys

if './' not in sys.path: sys.path.append('./')
from root.config.main import MPI, rAnk, mAster_rank

passed_nCSCG_RF2_tests = 0

start_time = MPI.Wtime()
if rAnk == mAster_rank: print(f"\n [nCSCG_RF2] tests start...\n")

from objects.nCSCG.rf2.base.__tests__.unittests.save_read import test_nCSCG_RF2_save_read


passed_nCSCG_RF2_tests += test_nCSCG_RF2_save_read()


if rAnk == mAster_rank:
    print("\n<{}> _nCSCG_RF2 tests passed; cost {:.3f} seconds.\n".format(
        passed_nCSCG_RF2_tests, MPI.Wtime()-start_time))
