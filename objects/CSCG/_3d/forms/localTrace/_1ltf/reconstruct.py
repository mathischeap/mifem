# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:55 PM
"""
from components.freeze.main import FrozenOnly


class _3dCSCG_1LocalTrace_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel=False, element_range=None):
        """

        Parameters
        ----------
        xi
        eta
        sigma
        ravel
        element_range : mesh element numbers

        Returns
        -------

        """
        raise NotImplementedError()

    def ___PrLT_make_reconstruction_matrix_on_grid___(self, xi, et, sg, element_range=None):
        """"""
        raise NotImplementedError()
