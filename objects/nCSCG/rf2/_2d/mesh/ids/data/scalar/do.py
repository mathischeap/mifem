# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 3:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_MeshIDS_Scalar_Do(FrozenOnly):
    """"""

    def __init__(self, scalar):
        """"""
        self._scalar_ = scalar
        self._freeze_self_()




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
