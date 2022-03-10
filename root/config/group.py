
from mpi4py import MPI
cOmm = MPI.COMM_WORLD
sIze: int = cOmm.Get_size()
rAnk: int = cOmm.Get_rank()
import numpy as np

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