# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 8:59 PM
"""
import sys

if './' not in sys.path: sys.path.append('../')

from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_MeshSpacePolynomialsBasisBase(FrozenOnly):
    """"""

    def __init__(self, mesh, xi_eta):
        """"""
        self._mesh_ = mesh
        self._xi_eta_ = xi_eta
        self._signature_ = mesh.signature

    @property
    def signature(self):
        return self._signature_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
