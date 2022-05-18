# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.do.find import mpRfT2_Mesh_Do_Find


class mpRfT2_Mesh_Do(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = mpRfT2_Mesh_Do_Find(self._mesh_)
        return self._find_










if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/do/main.py
    from objects.nCSCG.rf2._2d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3], 2, EDM='chaotic', show_info=True)
    mesh.do.unlock()

    if 0 in mesh.cscg.elements:
        c0 = mesh(0)
        c0.do.refine()

    mesh.do.update()