
from mpi4py import MPI
cOmm = MPI.COMM_WORLD
sIze: int = cOmm.Get_size()
rAnk: int = cOmm.Get_rank()
mAster_rank: int = 0 # you can, but you do not need to change this!
"""(int) the master core is under rank?"""
if sIze >= 2:
    sEcretary_rank: int = 1 # you can, but you do not need to change this!
    assert mAster_rank != sEcretary_rank, \
        f"when we have more than one core, master core should be different from secretary core."
else:
    sEcretary_rank: int = 0


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