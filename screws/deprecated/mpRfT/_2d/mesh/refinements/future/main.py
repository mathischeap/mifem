# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 8:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.refinements.future.preview import mpRfT2_Mesh_FutureRefinements_Preview
from objects.mpRfT._2d.mesh.refinements.future.h_adapt.main import mpRfT2_Mesh_FutureRefinements_hAdapt

class mpRfT2_Mesh_FutureRefinements(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._rfd_ = None
        self._preview_ = None
        self._h_adapt_to_ = mpRfT2_Mesh_FutureRefinements_hAdapt(mesh)
        self._freeze_self_()

    @property
    def rfd(self):
        return self._rfd_

    @property
    def preview(self):
        if self._preview_ is None:
            self._preview_ = mpRfT2_Mesh_FutureRefinements_Preview(self)
        return self._preview_
    
    @property
    def h_adapt_to(self):
        return self._h_adapt_to_


    





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/refinements/future/main.py
    pass
