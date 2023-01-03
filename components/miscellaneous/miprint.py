# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 5:32 PM
"""
from root.config.main import RANK, MASTER_RANK


def miprint(p, flush=True):

    if RANK == MASTER_RANK:
        print(p, flush=flush)
