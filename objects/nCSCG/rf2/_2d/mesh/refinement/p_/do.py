# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 2:34 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_Mesh_p_Refinement_Do(FrozenOnly):
    """"""

    def __init__(self, p):
        """"""
        self._p_ = p
        self._freeze_self_()


    def apply(self):
        self._p_.mesh.do.digest(self._p_)


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
