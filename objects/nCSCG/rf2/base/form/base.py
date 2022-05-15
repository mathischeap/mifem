# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/12 6:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class nCSCG_FormBase(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._k_ = None
        self.___is_wrapped_in_ADF___ = False # change this when we wrap it into an ADF

        self._freeze_self_()

    @property
    def mesh(self):
        return self._mesh_


    @property
    def k(self):
        return self._k_

    @property
    def ndim(self):
        """Return the dimensions."""
        return self.mesh.ndim




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
