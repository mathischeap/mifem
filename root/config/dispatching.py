
from mpi4py import MPI
cOmm = MPI.COMM_WORLD
sIze: int = cOmm.Get_size()
rAnk: int = cOmm.Get_rank()
mAster_rank: int = 0 # you can, but you do not need to change this!
sLave_ranks: list = [i for i in range(sIze)]
"""(list) The collections of ranks of all slaves (cores other than the master). So, 
the secretary core can be also a slave if ``sEcretary_rank != mAster_rank``."""
sLave_ranks.remove(mAster_rank)
import numpy as np




def dIspatching(originalTaskInputs: list, method):
    """
    Dispatch each input of ``originalTaskInputs`` from master to slaves.

    This is useful when the inputs and the results are not heavy. Results are not
    collected into the master.

    :param list originalTaskInputs: Only need data in master, be ``None`` in slaves.
    :param method: The method to take the input and do the job.
    """
    resultPOOL = []

    if rAnk == mAster_rank:
        taskInputs = originalTaskInputs.copy()
    else:
        taskInputs = None

    if sIze <= 2: # master does all jobs
        if rAnk == mAster_rank:
            while len(taskInputs) > 0:
                nextData = taskInputs.pop(0)
                resultPOOL.append(method(nextData))
        else:
            pass
    else:
        if rAnk == mAster_rank: # master does nothing, but distributing the tasks.
            CORE_STATUS = [True for _ in range(sIze)] # when True, the core is free.
            CORE_STATUS[mAster_rank] = False # Master core always not free.

            while len(taskInputs) > 0:

                if any(CORE_STATUS): # there are free cores.
                    core2bUsed = CORE_STATUS.index(True)
                    data2bSent = taskInputs.pop(0)
                    cOmm.send(data2bSent, dest=core2bUsed, tag=core2bUsed)
                    CORE_STATUS[core2bUsed] = False
                else:
                    newFreeCore = np.empty(1, dtype='i')
                    cOmm.Recv(newFreeCore, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
                    CORE_STATUS[newFreeCore[0]] = True

            for i in sLave_ranks:
                if CORE_STATUS[i]:
                    # it is free, so we do nothing
                    pass
                else:
                    # it is not free, then we know it is doing something, we have to wait for it to
                    # send back its rank, which means it is done with current job.
                    newFreeCore = np.empty(1, dtype='i')
                    cOmm.Recv(newFreeCore, source=i, tag=i)
                    assert newFreeCore[0] == i # make sure we recv something correct.
                    CORE_STATUS[newFreeCore[0]] = True

            # now all tasks are gone. Send message to all slaves to tell they stop waiting for more new tasks.
            for i in sLave_ranks:
                assert CORE_STATUS[i] # must be true
                cOmm.send('stop waiting for new tasks', dest=i, tag=i) # ask it stopping waiting.
        else:
            while 1:
                nextTaskData = cOmm.recv(source=mAster_rank, tag=rAnk)
                if nextTaskData == 'stop waiting for new tasks':
                    break # stop waiting
                else:
                    resultPOOL.append(method(nextTaskData))
                    cOmm.Send(np.array([rAnk], dtype='i'), dest=mAster_rank, tag=rAnk)

    cOmm.Barrier() # sync all cores.
    
    return resultPOOL
