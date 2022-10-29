# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/13 9:11 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT_Form_PartialCochain_Compile(FrozenOnly):
    """"""

    def __init__(self, pc):
        """"""
        self._f_ = pc._f_
        self._mesh_ = pc._mesh_
        self._pc_ = pc
        self._freeze_self_()

    @property
    def local(self):
        """The partial dofs are compiled as root-cell-wise local-dofs.

        Return a dict whose keys are the repr of the root-cells, and values are the local dofs.
        """
        dofs = self._pc_._pd_.compile.local
        CC = self._pc_._f_.cochain.__class__

        A = self._pc_._f_._cochain_

        self._pc_._f_._cochain_ = CC(self._pc_._f_)
        self._pc_._f_.discretization(target='boundary_condition')
        LC = self._pc_._f_.cochain.local
        cochains =  dict()
        for rc_rp in dofs:
            cochains[rc_rp] = LC[rc_rp][dofs[rc_rp]]

        self._pc_._f_._cochain_ = A
        return dofs, cochains







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/base/forms/base/partial/cochain/compile.py
    from __init__ import rfT2

    fc = rfT2.rf(100, mesh_pool='crazy')

    f1 = fc('1-f-o')

    import numpy as np
    def p(t, x, y): return np.sin(np.pi*x) * np.cos(np.pi*y) + t
    def q(t, x, y): return np.cos(np.pi*x) * np.sin(np.pi*y) + t

    v = fc('vector', (p, q))

    f1.BC.analytic_expression = v
    v.current_time = 0
    f1.BC.valid_boundaries = ['Upper', 'Left']
    fpc = f1.BC.partial_cochain

    D, C = fpc.compile.local
    print(D, C)


    t1 = fc('nst')
    s = fc('scalar', p)
    t1.BC.analytic_expression = s
    s.current_time = 0
    t1.BC.valid_boundaries = ['Upper', 'Left']
    tpc = t1.BC.partial_cochain

    D, C = tpc.compile.local
    print(D, C)