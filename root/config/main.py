# -*- coding: utf-8 -*-
"""
v3.1.0

In this script, we DO NOT use the structure of naming files and folders of the mifem library.

"""

# Version

vErsion = '3.0.1'


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
    seNtry_on = False # NEVER turn on this one: because we only monitor it through the master!


from root.config.tree import tRee
from root.config.dispatching import dIspatching
from root.config.chaining import cHaining
from root.config.group import gRoup_cores
from root.config.check_same_in_all_cores import cHeck_same_in_all_cores





if __name__ == '__main__':
    np = np
    tree = tRee
    d = dIspatching
    c = cHaining
    gc = gRoup_cores
    check = cHeck_same_in_all_cores