# -*- coding: utf-8 -*-
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE: int = COMM.Get_size()
RANK: int = COMM.Get_rank()
MASTER_RANK: int = 0  # you can, but you do not need to change this!
SLAVE_RANKS: list = [i for i in range(SIZE)]
"""(list) The collections of ranks of all slaves (cores other than the master). So, 
the secretary core can be also a slave if ``SECRETARY_RANK != MASTER_RANK``."""
SLAVE_RANKS.remove(MASTER_RANK)
import numpy as np


def DISPATCHING(originalTaskInputs: list, method):
    """
    Dispatch each input of ``originalTaskInputs`` from master to slaves.

    This is useful when the inputs and the results are not heavy. Results are not
    collected into the master.

    :param list originalTaskInputs: Only need data in master, be ``None`` in slaves.
    :param method: The method to take the input and do the job.
    """
    resultPOOL = []

    if RANK == MASTER_RANK:
        taskInputs = originalTaskInputs.copy()
    else:
        taskInputs = None

    if SIZE <= 2:  # master does all jobs
        if RANK == MASTER_RANK:
            while len(taskInputs) > 0:
                nextData = taskInputs.pop(0)
                resultPOOL.append(method(nextData))
        else:
            pass
    else:
        if RANK == MASTER_RANK:  # master does nothing, but distributing the tasks.
            CORE_STATUS = [True for _ in range(SIZE)]  # when True, the core is free.
            CORE_STATUS[MASTER_RANK] = False  # Master core always not free.

            while len(taskInputs) > 0:

                if any(CORE_STATUS):  # there are free cores.
                    core2bUsed = CORE_STATUS.index(True)
                    data2bSent = taskInputs.pop(0)
                    COMM.send(data2bSent, dest=core2bUsed, tag=core2bUsed)
                    CORE_STATUS[core2bUsed] = False
                else:
                    newFreeCore = np.empty(1, dtype='i')
                    COMM.Recv(newFreeCore, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
                    CORE_STATUS[newFreeCore[0]] = True

            for i in SLAVE_RANKS:
                if CORE_STATUS[i]:
                    # it is free, so we do nothing
                    pass
                else:
                    # it is not free, then we know it is doing something, we have to wait for it to
                    # send back its rank, which means it is done with current job.
                    newFreeCore = np.empty(1, dtype='i')
                    COMM.Recv(newFreeCore, source=i, tag=i)
                    assert newFreeCore[0] == i  # make sure we recv something correct.
                    CORE_STATUS[newFreeCore[0]] = True

            # now all tasks are gone. Send message to all slaves to tell they stop waiting for more new tasks.
            for i in SLAVE_RANKS:
                assert CORE_STATUS[i]  # must be true
                COMM.send('stop waiting for new tasks', dest=i, tag=i)  # ask it to stop waiting.
        else:
            while 1:
                nextTaskData = COMM.recv(source=MASTER_RANK, tag=RANK)
                if nextTaskData == 'stop waiting for new tasks':
                    break  # stop waiting
                else:
                    resultPOOL.append(method(nextTaskData))
                    COMM.Send(np.array([RANK], dtype='i'), dest=MASTER_RANK, tag=RANK)

    COMM.Barrier()  # sync all cores.
    
    return resultPOOL
