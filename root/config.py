# -*- coding: utf-8 -*-
"""
"""

# following value can be customized ...

caChe_factor = 2
"""We cache the intermediate data when there are ``caChe_factor`` repeated copies."""

assert isinstance(caChe_factor, int) and caChe_factor >= 1

saFe_mode = False
"""If sAfeMode is on, we will run many more checks."""

import numpy as np
from mpi4py import MPI
cOmm = MPI.COMM_WORLD

rAnk: int = cOmm.Get_rank()
"""(int) The rank I am currently in."""

sIze: int = cOmm.Get_size()
"""(int) How many cores in total we have?"""

mAster_rank: int = 0 # you can, but you do not need to change this!
"""(int) the master core is under rank?"""


assert sIze > 0, "I must have at least one core, right?"
if sIze >= 2:
    sEcretary_rank: int = 1 # you can, but you do not need to change this!
    assert mAster_rank != sEcretary_rank, \
        f"when we have more than one core, master core should be different from secretary core."
else:
    sEcretary_rank: int = 0


sLave_ranks: list = [i for i in range(sIze)]
"""(list) The collections of ranks of all slaves (cores other than the master). So, 
the secretary core can be also a slave if ``sEcretary_rank != mAster_rank``."""
sLave_ranks.remove(mAster_rank)

wOrker_ranks: list = [i for i in range(sIze)]
"""(list) The collections of ranks of all workers (cores other than the secretary). 
So, the master core can be also a worker if ``sEcretary_rank != mAster_rank``."""
wOrker_ranks.remove(sEcretary_rank)

assert 0 <= sEcretary_rank < sIze
assert 0 <= mAster_rank < sIze

# sentry setting, we can only turn on sentry in the master core!
if rAnk == mAster_rank:
    seNtry_on = False
    """If sEntryOn is on, we will monitor the scheme with Sentry."""
else:
    seNtry_on = False # NEVER turn on this one: cause we only monitor it through the master!



def tRee(factor=2):
    """We distribute our cores with a tree structure of branch ``factor``.

    :param int factor: (`default`: ``2``).
    :return: A generator.
    """
    if factor == 2:
        H = np.array([(i, i+1) if i+1<sIze else (i,-1) for i in range(0, sIze, 2)])
        yield ___parse_senders_and_receivers_factor2___(H)
        KEYS = H[:,0]
        LEN = np.size(KEYS)
        while LEN > 1:
            H = np.array([(KEYS[i], KEYS[i + 1]) if i + 1 < LEN else (KEYS[i], -1) for i in range(0, LEN, 2)])
            yield  ___parse_senders_and_receivers_factor2___(H)
            KEYS = H[:,0]
            LEN = np.size(KEYS)
    else:
        raise Exception(f"factor={factor} not coded.")

def ___parse_senders_and_receivers_factor2___(H):
    """"""
    for m,row in enumerate(H):
        if row[1] == -1:
            pass
        else:
            i, j = row
            i = int(i)
            j = int(j)
            tag = m
            if rAnk == i:
                recv_send = ('recv', {'source':j, 'tag':tag})
                return recv_send
            if rAnk == j:
                recv_send = ('send', {'dest':i, 'tag':tag})
                return recv_send


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

            # now all tasks are gone. Send message to slaves to tell they can stop waiting for more new tasks.
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


def cHaining(method, *args, **kwargs):
    """Let all cores do thing in a sequence; so next core will not start unless the previous core has
    sent him a message to start.

    :param method: method or function to be executed.
    :param args: Each core needs to have a copy.
    :param kwargs: Each core needs to have a copy.
    :return:
    """
    if sIze == 1:
        return method(*args, **kwargs)
    else:
        if rAnk == 0:
            R0 = method(*args, **kwargs)
            cOmm.send(1, dest=1, tag=0)
            return R0
        elif rAnk < sIze-1: # intermediate cores
            _ = cOmm.recv(source=rAnk-1, tag=rAnk-1)
            Ri = method(*args, **kwargs)
            cOmm.send(1, dest=rAnk+1, tag=rAnk)
            return Ri
        elif rAnk == sIze-1: # the last core
            _ = cOmm.recv(source=rAnk-1, tag=rAnk-1)
            return method(*args, **kwargs)
        else:
            raise Exception()

def gRoup_cores(member_num, group_num=None):
    """
    Group cores into groups each of ``member_num`` cores.

    :param member_num:
    :param group_num: When provide group_num, member_num does not matter at all.
    :return: For example, {'group number': 2, 'leaders': [0, 2], 0: [1], 2: [3]},
        which means we have 2 groups, and their leaders are core #0 and #2. The group
        lead by #0 has another member, core #1. Same for the group lead by core #2.
    """
    if group_num is None:
        if sIze == 1:
            group_num = 1
        else:
            if member_num > sIze:
                member_num = sIze
            elif member_num <= 1:
                member_num = 1
            else:
                member_num = int(member_num)
            group_num = int(np.ceil(sIze / member_num))
    else:
        pass
    assert 1 <= group_num <= sIze, f"group_num={group_num} wrong, only have {sIze} cores."
    num_cores = [sIze // group_num + (1 if x < sIze % group_num else 0) for x in range(group_num)]
    GROUPS = list()
    for i in range(group_num):
        GROUPS.append(range(sum(num_cores[0:i]), sum(num_cores[:(i + 1)])))
    assert max(GROUPS[-1]) == sIze-1, f"GROUPS={GROUPS} wrong."

    DG = dict()
    DG['group number'] = len(GROUPS)
    DG['leaders'] = list()

    for i, g in enumerate(GROUPS):
        assert len(g) > 0, f"group[{i}] is empty!"
        leader = min(g) # leader of group must be the core of smallest rank.
        DG[leader] = [i for i in g][1:]
        DG['leaders'].append(leader)
        if rAnk in g:
            DG['my leader'] = leader

    return DG

def cHeck_same_in_all_cores(*args):
    """
    check if some variables are the same in all cores.

    :param args:
    :return:
    """
    ARGS = cOmm.gather(args, root=sEcretary_rank)
    if rAnk == sEcretary_rank:
        ASSERTION = list()
        for i in range(sIze):
            ASSERTION.append(ARGS[i] == ARGS[0])
        ToF = all(ASSERTION)
        not_same_cores=list()
        for i in range(sIze):
            if not ASSERTION[i]:
                not_same_cores.append(i)
    else:
        ToF = None
        not_same_cores = None

    ToF, not_same_cores = cOmm.bcast([ToF,not_same_cores], root=sEcretary_rank)

    assert ToF, f"inputs in all cores must be same. These cores: {not_same_cores} have different inputs to core 0."
