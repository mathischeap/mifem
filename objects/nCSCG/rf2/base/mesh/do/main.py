# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/06 2:56 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

class nCSCG_MeshDoBase(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh

    def lock(self):
        """Can not refine or dulite the mesh!"""
        self._mesh_._locker_ = True

    def unlock(self):
        """Can refine or dulite the mesh!"""
        self._mesh_._locker_ = False







if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
