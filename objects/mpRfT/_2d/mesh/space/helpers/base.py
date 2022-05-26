# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/24/2022 10:53 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_Mesh_Space_Basis(FrozenOnly):
    """"""
    def __init__(self, mesh, coo_map):
        self._mesh_ = mesh
        if coo_map.distribution == 'uniform':
            assert coo_map._ndim_ == 1, f"evaluate_basis only accepts 1-d inputs."
            self._getitem_ = self.___Pr_getitem_uniform___
        elif coo_map.distribution == 'Gauss':
            self._getitem_ = self.___Pr_getitem_Gauss___
        else:
            raise NotImplementedError(f"{coo_map.__class__.__name__}")
        self._cm_ = coo_map
        self._cache_ = dict()

    def __getitem__(self, rp):
        raise NotImplementedError()

    def ___Pr_getitem_uniform___(self, rp):
        raise NotImplementedError()

    def ___Pr_getitem_Gauss___(self, rp):
        raise NotImplementedError()

    @property
    def ___Pr_rcMC_key___(self):
        """"""
        return self._cm_.___Pr_rcMC_key___

    def ___Pr_rcMC_nodes___(self, rp):
        """"""
        return self[rp][0]



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
