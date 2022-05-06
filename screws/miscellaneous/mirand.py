# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 3:42 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')


import random
from root.config.main import rAnk, mAster_rank, cOmm


def randint(a, b):
    """"""
    if rAnk == mAster_rank:
        r = random.randint(a, b)
    else:
        r = None
    return cOmm.bcast(r, root=mAster_rank)







if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
