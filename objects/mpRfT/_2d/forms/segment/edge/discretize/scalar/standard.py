# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

class mpRfT2_ESgF_Discretize_Standard_Scalar(FrozenOnly):
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

        F = F.___Pr_evaluate_func___()[0]







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/node/discretize/scalar/standard.py
    from __init__ import rfT2

    fc = rfT2.rf(100, N_range=(2,2))

    t0 = fc('nst', ndp=0, ntype='Lobatto')

    mesh = t0.mesh

    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar
    import numpy as np
    def p(t, x, y): return 1 / np.exp(np.abs(x**2 + y**2 - 2)) + t
    s = mpRfT2_Scalar(mesh, p)

    t0.TW.func = s
    s.current_time = 0
    t0.discretize()



    t0.reconstruct()