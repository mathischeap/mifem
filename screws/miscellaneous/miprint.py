# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 5:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from root.config.main import rAnk, mAster_rank

def miprint(p, flush=True):

    if rAnk == mAster_rank:
        print(p, flush=flush)


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
