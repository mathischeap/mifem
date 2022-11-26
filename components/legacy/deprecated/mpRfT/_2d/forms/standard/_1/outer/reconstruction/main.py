# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 5:08 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
import numpy as np
from objects.mpRfT._2d.forms.standard._1.outer.reconstruction._matrix import mpRfT2_So1F_Reconstruction_tMatrix
from objects.mpRfT._2d.forms.standard._1.outer.reconstruction.matrix import mpRfT2_So1F_Reconstruction_Matrix


class mpRfT2_So1F_Reconstruction(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, coo_map, rp=None, ravel=False, value_only=False):
        """

        Parameters
        ----------
        coo_map
        ravel
        rp
        value_only

        Returns
        -------

        """
        assert coo_map.___Pr_is_mpRfT2_mesh_coo_map___, f"I need a mpRfT2 coo_map object."
        f = self._f_
        mesh = f.mesh

        # ---- parse indices -----------------------------------------------------------
        if rp is None:
            rps = mesh.rcfc
        else:
            raise NotImplementedError(f"rp={rp} is not implemented.")
        #================================================================================

        if rps is mesh.rcfc:
            full = True
        else:
            full = False

        Basis = mesh.space.do.evaluate_basis(self._f_, coo_map)

        rcC_invJ = mesh.rcMC.inverse_Jacobian_matrix(Basis)

        if value_only:
            value = dict()

            for rp in rps:
                xi_et, basis, _2d_shape = Basis[rp]
                cell= mesh[rp]
                invJ = rcC_invJ[rp]

                LOCALS = f.cochain.___Pr_divide_local___(rp)
                u = np.einsum('ij, i -> j', basis[0], LOCALS[0], optimize='greedy')
                v = np.einsum('ij, i -> j', basis[1], LOCALS[1], optimize='greedy')

                vi = [None, None]
                typeWr2Metric = cell.type_wrt_metric.mark
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    vi[0] = u*invJ[1][1]
                    vi[1] = v*invJ[0][0]
                else:
                    vi[0] = + u*invJ[1][1] - v*invJ[0][1]
                    vi[1] = - u*invJ[1][0] + v*invJ[0][0]

                if ravel:
                    value[rp] = vi
                else:
                    value[rp] = [vi[0].reshape(_2d_shape, order='F'),
                                 vi[1].reshape(_2d_shape, order='F')]

            value = mesh.rcWds.vector(value, 2, coo_map.distribution, full)

            return value


        else:
            xy = dict()
            value = dict()

            for rp in rps:
                xi_et, basis, _2d_shape = Basis[rp]
                cell= mesh[rp]
                invJ = rcC_invJ[rp]
                xyi = cell.coordinate_transformation.mapping(*xi_et)

                LOCALS = f.cochain.___Pr_divide_local___(rp)
                u = np.einsum('ij, i -> j', basis[0], LOCALS[0], optimize='greedy')
                v = np.einsum('ij, i -> j', basis[1], LOCALS[1], optimize='greedy')

                vi = [None, None]
                typeWr2Metric = cell.type_wrt_metric.mark
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    vi[0] = u*invJ[1][1]
                    vi[1] = v*invJ[0][0]
                else:
                    vi[0] = + u*invJ[1][1] - v*invJ[0][1]
                    vi[1] = - u*invJ[1][0] + v*invJ[0][0]

                if ravel:
                    xy[rp] = xyi
                    value[rp] = vi
                else:
                    xy[rp] = [xyi[0].reshape(_2d_shape, order='F'),
                              xyi[1].reshape(_2d_shape, order='F')]
                    value[rp] = [vi[0].reshape(_2d_shape, order='F'),
                                 vi[1].reshape(_2d_shape, order='F')]

            xy = mesh.rcWds.vector(xy, 2, coo_map.distribution, full)
            value = mesh.rcWds.vector(value, 2, coo_map.distribution, full)

            return xy, value

    @property
    def ___Pr_test_matrix___(self):
        return mpRfT2_So1F_Reconstruction_tMatrix(self._f_)

    @property
    def matrix(self):
        return mpRfT2_So1F_Reconstruction_Matrix(self._f_)





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/outer/reconstruct.py
    from __init__ import rfT2
    fc = rfT2.rf(100, N_range=(3,4))
    # from objects.mpRfT._2d.master import MeshGenerator, FormCaller
    # mesh = MeshGenerator('crazy', c=0.1)([20, 20], 3, show_info=True)
    # fc = FormCaller(mesh)

    f = fc('1-f-o')

    def p(t, x, y): return np.sin(np.pi*x) * np.cos(np.pi*y) + t
    def q(t, x, y): return np.cos(np.pi*x) * np.sin(np.pi*y) + t
    v = fc('vector', (p, q))

    f.TW.func = v
    v.current_time = 0

    f.discretize()

    coo = f.mesh.coo_map.uniform(20, ndim=1)

    f.visualize(show_mesh=True)