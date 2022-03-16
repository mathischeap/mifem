# -*- coding: utf-8 -*-
"""
Here we test our iterators.

$ mpiexec -n 2 python tests\iterator_tests.py

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from time import sleep

from tools.iterators.simple import SimpleIterator


def ___TEST_SOLVER___(tk, tk1):
    """
    Parameters
    ----------
    tk :
    tk1 :

    Returns
    -------
    exit_code : The standard exit code.
    shut_down : If it is ``True``, the outer iterator will shut down immediately.
    message : The solver message.
    output1 :
    """
    sleep(1)
    assert tk1 > tk
    if rAnk == mAster_rank:
        ot1 = (tk+tk1)/2
    else:
        ot1 = None
    ot1 = cOmm.bcast(ot1, root=mAster_rank)
    message = ['Inner solver: [gmres] A shape=(92000, 92000), Solver: GMRES, routine: MPI1, restart=150, convergence info=100, residual=1e-06, number of iterations=100, time cost=390.98 seconds.',
               "Outer solver: [icpsNS-OS1] Total system shape (226400, 226400), Reduced system shape (103200, 103200) @sparsity~=99.86% with ~ 14967722 non-zeros. Solver costs 458.77s among which 418.69s is for solving the reduced system. [gmres0] on the reduced system: restart=150, convergence info=100, residual=1e-06, number of iterations=100"]
    return 1, 0, message, ot1


SI = SimpleIterator(t0=0, dt=0.1, max_steps=100000,
                    monitor_factor=1,
                    RDF_filename='RDF_filename',
                    name='test-iterator')


SI(___TEST_SOLVER___, [1,])

SI.run()

# if rAnk == mAster_rank: os.remove('RDF_filename.csv')


# S2 = SimpleIterator.read('RDF_filename')
#
#
# if rAnk == mAster_rank:
#     print(S2._solver_source_code_)
#     print(S2._solver_)
#     print(S2._solver_dir_)
#
#
# if rAnk == mAster_rank: os.remove('RDF_filename.mitr')