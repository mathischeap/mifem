# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/01 3:13 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_Mesh_sgMC_Base(FrozenOnly):
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
        self._key_ = coo.___Pr_sgMC_key___
        self._nodes_ = coo.___Pr_sgMC_nodes___

        return self

    def __getitem__(self, sg):
        raise NotImplementedError()




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
