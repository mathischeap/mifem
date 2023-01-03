# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 3:42 PM
"""
import random
from root.config.main import RANK, MASTER_RANK, COMM


def randint(a, b):
    """"""
    if RANK == MASTER_RANK:
        r = random.randint(a, b)
    else:
        r = None
    return COMM.bcast(r, root=MASTER_RANK)


def sample(population, k):
    if RANK == MASTER_RANK:
        r = random.sample(population, k)
    else:
        r = None
    return COMM.bcast(r, root=MASTER_RANK)
