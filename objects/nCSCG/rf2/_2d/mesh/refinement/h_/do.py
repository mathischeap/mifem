# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 2:34 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_Mesh_h_Refinement_Do(FrozenOnly):
    """"""

    def __init__(self, h):
        """"""
        self._h_ = h
        self._freeze_self_()


    def apply(self):
        self._h_.mesh.do.digest(self._h_)




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
