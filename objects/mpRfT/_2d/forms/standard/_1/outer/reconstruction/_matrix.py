# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/22 5:50 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np
from scipy.sparse import csr_matrix


class mpRfT2_So1F_Reconstruction_tMatrix(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._coo_map_ = None
        self._basis_ = None
        self._freeze_self_()

    def __Initialize__(self, coo_map):
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

    def __call__(self, rc_rp):
        """Return the reconstruction matrix in the root-cell: rc_rp."""
        basis = self._basis_[rc_rp]

        cell = self._f_.mesh[rc_rp]
        mark = cell.type_wrt_metric.mark

        COO = self._coo_map_[rc_rp]
        u_ouv = cell.coordinate_transformation.unit_outward_norm_vector()

        VV = list()
        for edge in basis:
            XI_ETA = basis[edge]
            for i, xi_eta in enumerate(XI_ETA):
                xi, eta = xi_eta[0]
                invJ = cell.coordinate_transformation.inverse_Jacobian_matrix(xi, eta)
                J = cell.coordinate_transformation.Jacobian(xi, eta)

                u, v = xi_eta[1]

                _U = np.vstack([+ u * invJ[1][1], - v * invJ[0][1]])
                _V = np.vstack([- u * invJ[1][0], + v * invJ[0][0]])

                if isinstance(mark, str) and mark[:4] == 'Orth':
                    _ = u_ouv[edge]
                    V = _U * _[0] + _V * _[1]
                    weights, intervals = COO[edge][i][2:4]
                    LEN = weights.shape[0]
                    VS = V.shape
                    _ = int(VS[1]/LEN)
                    V = V.reshape((VS[0], _, _), order='F')
                    V = np.einsum('ijk, j, k  -> ki', V, weights, intervals*0.5*(J**0.5)[0])


                else:
                    raise NotImplementedError()

                VV.append(V)

        VV = np.vstack(VV)
        VV[abs(VV)<1e-10] = 0

        return csr_matrix(VV)





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
