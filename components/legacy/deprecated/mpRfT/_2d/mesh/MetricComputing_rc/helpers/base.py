# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/24/2022 9:43 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_Mesh_rcMC_Base(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._cache_ = dict()

        self._coo_ = None
        self._key_ = None
        self._nodes_ = None

        self._freeze_self_()

    def __call__(self, coo):
        """"""
        assert coo._mesh_ is self._mesh_
        self._coo_ = coo
        self._key_ = coo.___Pr_rcMC_key___
        self._nodes_ = coo.___Pr_rcMC_nodes___
        return self

    def __getitem__(self, rp):
        raise NotImplementedError()




if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/rcMetricComputing/helpers/base.py
    pass
