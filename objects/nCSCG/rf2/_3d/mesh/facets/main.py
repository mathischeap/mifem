# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12:23 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._3d.mesh.facets.do import FacetsDo

class _3nCSCG_MeshFacets(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._do_ = None
        self._freeze_self_()

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = FacetsDo(self)
        return self._do_



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
