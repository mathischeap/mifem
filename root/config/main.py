# -*- coding: utf-8 -*-
"""
v3.1.0

In this script, we DO NOT use the structure of naming files and folders of the mifem library.

"""

# Version

VERSION = '3.0.1'
# following value can be customized ...

CACHE_FACTOR = 2
"""We cache the intermediate data when there are ``CACHE_FACTOR`` repeated copies."""

assert isinstance(CACHE_FACTOR, int) and CACHE_FACTOR >= 1

SAFE_MODE = False
"""If SAFE_MODE is on, we will do many more checks."""

import numpy as np
from mpi4py import MPI
COMM = MPI.COMM_WORLD

RANK: int = COMM.Get_rank()
"""(int) The rank I am currently in."""

SIZE: int = COMM.Get_size()
"""(int) How many cores in total we have?"""

MASTER_RANK: int = 0  # you can, but you do not need to change this!
"""(int) the master core is under rank?"""

if RANK == MASTER_RANK:
    ASSEMBLE_COST = {
        'accumulated': 0,  # we have used this much time (seconds) on assembling.
        'recent': list(),  # we additionally have spent this much time on assembling.
    }
else:
    ASSEMBLE_COST = None

assert SIZE > 0, "I must have at least one core, right?"
if SIZE >= 2:
    SECRETARY_RANK: int = 1  # you can, but you do not need to change this!
    assert MASTER_RANK != SECRETARY_RANK, \
        f"when we have more than one core, master core should be different from secretary core."
else:
    SECRETARY_RANK: int = 0

SLAVE_RANKS: list = [i for i in range(SIZE)]
"""(list) The collections of ranks of all slaves (cores other than the master). So, 
the secretary core can be also a slave if ``SECRETARY_RANK != MASTER_RANK``."""
SLAVE_RANKS.remove(MASTER_RANK)

WORKER_RANKS: list = [i for i in range(SIZE)]
"""(list) The collections of ranks of all workers (cores other than the secretary). 
So, the master core can be also a worker if ``SECRETARY_RANK != MASTER_RANK``."""
WORKER_RANKS.remove(SECRETARY_RANK)

assert 0 <= SECRETARY_RANK < SIZE
assert 0 <= MASTER_RANK < SIZE

# sentry setting, we can only turn on sentry in the master core!
if RANK == MASTER_RANK:
    SENTRY_ON = False
    """If SENTRY_ON is True, we will monitor the scheme with Sentry."""
else:
    SENTRY_ON = False  # NEVER turn on this one because we only monitor it through the master!


from root.config.tree import TREE
from root.config.dispatching import DISPATCHING
from root.config.chaining import CHAINING
from root.config.group import GROUP_CORES
from root.config.check_same_in_all_cores import CHECK_SAME_IN_ALL_CORES


if __name__ == '__main__':
    np = np
    tree = TREE
    d = DISPATCHING
    c = CHAINING
    gc = GROUP_CORES
    check = CHECK_SAME_IN_ALL_CORES
