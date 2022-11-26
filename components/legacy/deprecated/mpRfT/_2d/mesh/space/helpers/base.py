# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/24/2022 10:53 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_Mesh_Space_Basis(FrozenOnly):
    """"""
    def __init__(self, mesh, coo_map):
        self._mesh_ = mesh
        if coo_map.distribution == 'uniform':
            assert coo_map._ndim_ == 1, f"evaluate_basis only accepts 1-d inputs."
            self._getitem_ = self.___Pr_getitem_uniform___
        elif coo_map.distribution == 'Gauss':
            self._getitem_ = self.___Pr_getitem_Gauss___
        elif coo_map.distribution == 'Lobatto':
            self._getitem_ = self.___Pr_getitem_Lobatto___
        elif coo_map.distribution == 'PSI':
            self._getitem_ = self.___Pr_getitem_PSI___
        elif coo_map.distribution == 'SegInt':
            self._getitem_ = self.___Pr_getitem_SegInt___

        else:
            raise NotImplementedError(f"{self.__class__.__name__} cannot handle {coo_map.__class__.__name__}")
        self._cm_ = coo_map
        self._cache_ = dict()

    def __getitem__(self, rp):
        raise NotImplementedError()

    def ___Pr_getitem_uniform___(self, rp):
        raise NotImplementedError()

    def ___Pr_getitem_Gauss___(self, rp):
        raise NotImplementedError()

    def ___Pr_getitem_Lobatto___(self, rp):
        raise NotImplementedError()

    def ___Pr_getitem_PSI___(self, rp):
        raise NotImplementedError()

    def ___Pr_getitem_SegInt___(self, rp):
        raise NotImplementedError()

    @property
    def ___Pr_rcMC_key___(self):
        """"""
        raise NotImplementedError()

    def ___Pr_rcMC_nodes___(self, rp):
        """"""
        raise NotImplementedError()

    @property
    def ___Pr_sgMC_key___(self):
        """"""
        raise NotImplementedError()

    def ___Pr_sgMC_nodes___(self, rp):
        """"""
        raise NotImplementedError()






if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
