# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 12:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from screws.freeze.base import FrozenOnly
from screws.quadrature import Quadrature


class mpRfT2_So1F_Discretize_Standard_Vector(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self.___DISCRETIZE_STANDARD_CACHE___ = dict()
        self._freeze_self_()

    def ___PRIVATE_discretize_preparation___(self, N):
        """

        Parameters
        ----------
        N

        Returns
        -------

        """
        if N in self.___DISCRETIZE_STANDARD_CACHE___:
            return self.___DISCRETIZE_STANDARD_CACHE___[N]

        self.___DISCRETIZE_STANDARD_CACHE___[N] = list()

        p = [N + 2, N + 2]

        quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
        p_x, p_y = p

        space = self._f_.mesh.space[N]

        edges_size = [space.nodes[i][1:] - space.nodes[i][:-1] for i in range(2)]

        cell_nodes = [(0.5 * (edges_size[i][np.newaxis, :]) * (quad_nodes[i][:, np.newaxis] + 1)
            + space.nodes[i][:-1]).ravel('F') for i in range(2)]


        quad_xi = np.tile(cell_nodes[0], N + 1).reshape(
            (p_x + 1, N * (N + 1)), order='F')
        quad_eta = np.repeat(space.nodes[1][np.newaxis, :], N, axis=0).ravel('F')
        quad_eta = quad_eta[np.newaxis, :].repeat(p_x + 1, axis=0)
        ES = np.tile(edges_size[0], N + 1)
        self.___DISCRETIZE_STANDARD_CACHE___[N].append([quad_xi, quad_eta, ES, quad_weights])

        quad_xi = np.tile(space.nodes[0], N)[np.newaxis, :].repeat(p_y + 1, axis=0)
        quad_eta = np.repeat(cell_nodes[1].reshape(
            (p_y + 1, N), order='F'), N + 1, axis=1)
        ES = np.repeat(edges_size[1], N + 1)
        self.___DISCRETIZE_STANDARD_CACHE___[N].append([quad_xi, quad_eta, ES, quad_weights])

        return self.___DISCRETIZE_STANDARD_CACHE___[N]


    def __call__(self, target):
        """

        Parameters
        ----------
        target : str
            {'func',}

        Returns
        -------

        """
        # --- target --------------------------------------------------------
        if target == 'func':
            FUNC = self._f_.TW.func
        else:
            raise NotImplementedError(f"I cannot deal with target = {target}.")

        FUNC = FUNC.___Pr_evaluate_func___()

        local_dx = dict()
        local_dy = dict()
        JXC, JYC = dict(), dict() # local caches

        rcfc = self._f_.mesh.rcfc
        for rp in rcfc:
            cell = rcfc[rp]
            N = cell.N

            DD = self.___PRIVATE_discretize_preparation___(N)
            xi_x, eta_x, edge_size_d_xi, quad_weights = DD[0]
            xi_y, eta_y, edge_size_d_et, quad_weights = DD[1]
            edge_size = (edge_size_d_xi, edge_size_d_et)

            typeWr2Metric = cell.___Pr_metric_N_key___

            xy = cell.coordinate_transformation.mapping(xi_x, eta_x)

            if typeWr2Metric in JXC:
                J = JXC[typeWr2Metric]
            else:
                J = cell.coordinate_transformation.Jacobian_matrix(xi_x, eta_x)

                if not typeWr2Metric.isnumeric():
                    JXC[typeWr2Metric] = J

            if typeWr2Metric[:4] == 'Orth':
                v = FUNC[1](*xy)
                local_dx[rp] = np.einsum(
                    'jk, j, k -> k', J[0][0]*v, quad_weights[0], edge_size[0] * 0.5, optimize='greedy'
                )

            else:
                J = (J[0][0], J[1][0])
                u = FUNC[0](*xy)
                v = FUNC[1](*xy)
                local_dx[rp] = np.einsum(
                    'jk, j, k -> k', J[0]*v - J[1]*u, quad_weights[0], edge_size[0] * 0.5, optimize='greedy'
                )

            xy = cell.coordinate_transformation.mapping(xi_y, eta_y)
            if typeWr2Metric in JYC:
                J = JYC[typeWr2Metric]
            else:
                J = cell.coordinate_transformation.Jacobian_matrix(xi_y, eta_y)
                if not typeWr2Metric.isnumeric():
                    JYC[typeWr2Metric] = J

            if typeWr2Metric[:4] == 'Orth':
                u = FUNC[0](*xy)
                local_dy[rp] = np.einsum(
                    'jk, j, k -> k', J[1][1]*u, quad_weights[1], edge_size[1]*0.5, optimize='greedy'
                )
            else:
                J = (J[0][1], J[1][1])
                u = FUNC[0](*xy)
                v = FUNC[1](*xy)
                local_dy[rp] = np.einsum(
                    'jk, j, k -> k', - J[0]*v + J[1]*u, quad_weights[1], edge_size[1]*0.5, optimize='greedy'
                )

        del JXC, JYC

        cochainLocal = dict()
        for rp in rcfc:
            cochainLocal[rp] = np.hstack((local_dy[rp], local_dx[rp]))

        return cochainLocal









if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/outer/discretize/vector/standard.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    f = fc('1-f-o')

    def p(t, x, y): return np.sin(np.pi*x) * np.cos(np.pi*y) + t
    def q(t, x, y): return np.cos(np.pi*x) * np.sin(np.pi*y) + t
    v = fc('vector', (p, q))

    f.TW.func = v
    v.current_time = 0

    f.discretize()
