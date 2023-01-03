# -*- coding: utf-8 -*-
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE: int = COMM.Get_size()
RANK: int = COMM.Get_rank()
MASTER_RANK: int = 0  # you can, but you do not need to change this!
"""(int) the master core is under rank?"""
if SIZE >= 2:
    SECRETARY_RANK: int = 1  # you can, but you do not need to change this!
    assert MASTER_RANK != SECRETARY_RANK, \
        f"when we have more than one core, master core should be different from secretary core."
else:
    SECRETARY_RANK: int = 0


def CHECK_SAME_IN_ALL_CORES(*args):
    """
    check if some variables are the same in all cores.

    :param args:
    :return:
    """
    ARGS = COMM.gather(args, root=SECRETARY_RANK)
    if RANK == SECRETARY_RANK:
        ASSERTION = list()
        for i in range(SIZE):
            ASSERTION.append(ARGS[i] == ARGS[0])
        ToF = all(ASSERTION)
        not_same_cores = list()
        for i in range(SIZE):
            if not ASSERTION[i]:
                not_same_cores.append(i)
    else:
        ToF = None
        not_same_cores = None

    ToF, not_same_cores = COMM.bcast([ToF, not_same_cores], root=SECRETARY_RANK)

    assert ToF, f"inputs in all cores must be same. These cores: {not_same_cores} have different inputs to core 0."
