# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/25/2022 9:35 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from screws.freeze.base import FrozenOnly


class mpRfT2_ESgF_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._freeze_self_()

    def __call__(self, coo_map, seg=None, value_only=False):
        """

        Parameters
        ----------
        coo_map
        seg :
            In which segments we are going to reconstruct the standard form.

        Returns
        -------

        """
        assert coo_map.___Pr_is_mpRfT2_mesh_coo_map___, f"I need a mpRfT2 coo_map object."
        t = self._t_
        mesh = t.mesh

        # ---- parse indices -----------------------------------------------------------
        if seg is None:
            segments = mesh.segments
        else:
            raise NotImplementedError(f"seg={seg} is not implemented.")
        #================================================================================

        Basis = mesh.space.do.evaluate_basis(self._t_, coo_map)






if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/node/reconstruct/main.py

    from __init__ import rfT2

    fc = rfT2.rf(100, N_range=(2,2))

    t0 = fc('nst', ndp=0, ntype='Lobatto')

    mesh = t0.mesh

    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar
    def p(t, x, y): return np.sin(np.pi*x) * np.sin(np.pi*y) + t
    s = mpRfT2_Scalar(mesh, p)

    t0.TW.func = s
    s.current_time = 0
    t0.discretize()

    coo = mesh.coo_map.uniform(20, ndim=1)

    xy, v = t0.reconstruct(coo)

    v.visualize(xy)