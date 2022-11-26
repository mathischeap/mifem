# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/13 9:11 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT_Form_PartialCochain_Include(FrozenOnly):
    """"""

    def __init__(self, pc):
        """"""
        self._f_ = pc._f_
        self._mesh_ = pc._mesh_
        self._pc_ = pc

        self._freeze_self_()

    def boundaries(self, boundary_names):
        """We will add keys to self._pd_.keys. Those added keys represents the surfaces of root-cells
        on the boundary names.
        """
        self._pc_._pd_.include.boundaries(boundary_names)




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/base/forms/base/partial/cochain/include.py
    from __init__ import rfT2

    fc = rfT2.rf(100, mesh_pool='crazy')

    f1 = fc('1-f-o')
    f1.BC.valid_boundaries = ['Upper', 'Left']
    fpc = f1.BC.partial_cochain

    import numpy as np
    def p(t, x, y): return np.sin(np.pi*x) * np.cos(np.pi*y) + t
    def q(t, x, y): return np.cos(np.pi*x) * np.sin(np.pi*y) + t

    v = fc('vector', (p, q))

    f1.TW.func = v
    v.current_time = 0

    print(fpc._pd_._indicators_)
