# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 5:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from root.config.main import RANK, MASTER_RANK

def miprint(p, flush=True):

    if RANK == MASTER_RANK:
        print(p, flush=flush)


if __name__ == "__main__":
    # mpiexec -n 4 python screws/miscellaneous/miprint.py
    pass
