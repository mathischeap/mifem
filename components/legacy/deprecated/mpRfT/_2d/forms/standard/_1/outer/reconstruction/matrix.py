# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 6/28/2022 10:11 AM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
import numpy as np
from scipy.sparse import csr_matrix


class mpRfT2_So1F_Reconstruction_Matrix(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._coo_map_ = None
        self._basis_ = None
        self._freeze_self_()

    def __call__(self, coo_map):
        """Also for the initialization."""
        assert coo_map.___Pr_is_mpRfT2_mesh_coo_map___, f"I need a mpRfT2 coo_map object."
        assert coo_map._mesh_ is self._f_.mesh
        self._coo_map_ = coo_map
        mesh = self._f_.mesh
        self._basis_ = mesh.space.do.evaluate_basis(self._f_, coo_map)
        return self

    def __iter__(self):
        for rc_rp in self._f_.mesh.rcfc:
            yield rc_rp

    def __getitem__(self, rc_rp):
        """Return the reconstruction matrix in the root-cell: rc_rp."""
        basis = self._basis_[rc_rp]

        cell = self._f_.mesh[rc_rp]
        mark = cell.type_wrt_metric.mark

        VD = dict()
        for edge in basis:
            XI_ETA = basis[edge]
            VE = list()
            for i, xi_eta in enumerate(XI_ETA):
                xi, eta = xi_eta[0]
                invJ = cell.coordinate_transformation.inverse_Jacobian_matrix(xi, eta)

                u, v = xi_eta[1]

                _U = np.vstack([+ u * invJ[1][1], - v * invJ[0][1]]).T
                _V = np.vstack([- u * invJ[1][0], + v * invJ[0][0]]).T

                if isinstance(mark, str) and mark[:4] == 'Orth':
                    if edge   == 'U': V = - _U
                    elif edge == 'D': V = + _U
                    elif edge == 'L': V = - _V
                    else:             V = + _V
                else:
                    raise NotImplementedError()

                VE.append(V)

            VD[edge] = VE

        return VD





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/outer/reconstruction/matrix.py
    from __init__ import rfT2
    fc = rfT2.rf(100, N_range=(3,4))
    # from objects.mpRfT._2d.master import MeshGenerator, FormCaller
    # mesh = MeshGenerator('crazy', c=0.1)([20, 20], 3, show_info=True)
    # fc = FormCaller(mesh)

    f = fc('1-f-o')

    def p(t, x, y): return np.sin(np.pi*x) * np.cos(np.pi*y) + t
    def q(t, x, y): return np.cos(np.pi*x) * np.sin(np.pi*y) + t
    v = fc('vector', (p, q))

    f.analytic_expression = v
    v.current_time = 0

    f.discretization()
    f.visualization()
