# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly

class mpRfT2_S0F_Discretize_Standard_Scalar(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
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
        mesh = self._f_.mesh

        if target == 'analytic_expression':
            F = self._f_.analytic_expression
        else:
            raise NotImplementedError()

        F = F.___Pr_evaluate_func___()[0]

        rcw_LC = dict() # root-cell-wise Local cochain

        for rp in mesh.rcfc: # go through all local root cells.
            cell = mesh[rp] # local root cell
            xi_et = cell.coordinate_transformation.mapping(*cell.space.GoN_ravel)
            rcw_LC[rp] = F(*xi_et)

        return rcw_LC










if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
