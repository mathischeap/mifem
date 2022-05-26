# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

class mpRfT2_T1F_Discretize_Standard_Scalar(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._freeze_self_()

    def __call__(self, target):
        """

        Parameters
        ----------
        target : str
            {'func',}

        Returns
        -------

        """
        mesh = self._t_.mesh

        if target == 'func':
            F = self._t_.TW.func
        else:
            raise NotImplementedError()

        # F = F.___Pr_evaluate_func___()[0]
        #
        # rcw_LC = dict() # root-cell-wise Local cochain
        #
        # for rp in mesh.rcfc: # go through all local root cells.
        #     cell = mesh[rp] # local root cell
        #     xi_et = cell.coordinate_transformation.mapping(*cell.space.GoN_ravel)
        #     rcw_LC[rp] = F(*xi_et)
        #
        # return rcw_LC










if __name__ == "__main__":
    # mpiexec -n 4 python 
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t1 = fc('1-t')
