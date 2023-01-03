# -*- coding: utf-8 -*-
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE: int = COMM.Get_size()
RANK: int = COMM.Get_rank()
import numpy as np


def GROUP_CORES(member_num, group_num=None):
    """
    Group cores into groups each of ``member_num`` cores.

    :param member_num:
    :param group_num: When provide group_num, member_num does not matter at all.
    :return: For example, {'group number': 2, 'leaders': [0, 2], 0: [1], 2: [3]},
        which means we have 2 groups, and their leaders are core #0 and #2. The group
        lead by #0 has another member, core #1. Same for the group lead by core #2.
    """
    if group_num is None:
        if SIZE == 1:
            group_num = 1
        else:
            if member_num > SIZE:
                member_num = SIZE
            elif member_num <= 1:
                member_num = 1
            else:
                member_num = int(member_num)
            group_num = int(np.ceil(SIZE / member_num))
    else:
        pass
    assert 1 <= group_num <= SIZE, f"group_num={group_num} wrong, only have {SIZE} cores."
    num_cores = [SIZE // group_num + (1 if x < SIZE % group_num else 0) for x in range(group_num)]
    GROUPS = list()
    for i in range(group_num):
        GROUPS.append(range(sum(num_cores[0:i]), sum(num_cores[:(i + 1)])))
    assert max(GROUPS[-1]) == SIZE - 1, f"GROUPS={GROUPS} wrong."

    DG = dict()
    DG['group number'] = len(GROUPS)
    DG['leaders'] = list()

    for i, g in enumerate(GROUPS):
        assert len(g) > 0, f"group[{i}] is empty!"
        leader = min(g)  # leader of group must be the core of the smallest rank.
        DG[leader] = [i for i in g][1:]
        DG['leaders'].append(leader)
        if RANK in g:
            DG['my leader'] = leader

    return DG
