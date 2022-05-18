# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2.base.form.standard.base.cochain.local.RCW_full import \
    nCSCG_RF2__RCW_Full__LocalCochain

class _2nCSCG_RF2_S0F_Discretize_Standard_Scalar(FrozenOnly):
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

        if target == 'func':
            F = self._f_.TW.func
        else:
            raise NotImplementedError()

        F = F.___Pr_evaluate_func___()[0]

        ScW_LC = dict() # Sub-cell-Wise Local cochain
        for indices in mesh: # go through all local root cells.
            cell = mesh(indices) # local root cell
            xi_et = cell.coordinate_transformation.mapping(*cell.space.GoN_ravel)
            ScW_LC[repr(cell)] = F(*xi_et)

        return nCSCG_RF2__RCW_Full__LocalCochain(self._f_, ScW_LC)










if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
