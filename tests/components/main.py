# -*- coding: utf-8 -*-
"""
To run all test with given number of threads, do

mpiexec -n 6 python tests/components/main.py

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/21/2022 3:55 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *

passed_components_tests = 0

if RANK == MASTER_RANK: print(f"\n [components] tests start...\n")

t_global_start = MPI.Wtime()


from tests.components.ft2dw import test_functions_time_plus_2d_wrappers_AKA_ft2dw
from tests.components.ft3dw import test_functions_time_plus_3d_wrappers_AKA_ft3dw
from tests.components.assemblers import test_vector_assembler



passed_components_tests += test_functions_time_plus_2d_wrappers_AKA_ft2dw()
passed_components_tests += test_functions_time_plus_3d_wrappers_AKA_ft3dw()
passed_components_tests += test_vector_assembler()



COMM.barrier()
if RANK == MASTER_RANK:
    print("\n<{}> [components] tests passed; cost {:.3f} seconds.\n".format(
        passed_components_tests, MPI.Wtime()-t_global_start))