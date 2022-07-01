# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np
from screws.freeze.base import FrozenOnly
from screws.quadrature import Quadrature

class mpRfT2_S2F_Discretize_Standard_Scalar(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self.___D_cache___ = dict()
        self._freeze_self_()


    def ___Pr_prepare_discretize_data___(self, N):
        """"""
        if N in self.___D_cache___:
            return self.___D_cache___[N]

        f = self._f_
        space = f.mesh.space[N]

        p = [N + 2, N + 2]
        quad_nodes, quad_weights = Quadrature(p).quad

        magic_factor = 0.25
        xi = np.zeros((N**2, p[0] + 1, p[1] + 1))
        et = np.zeros((N**2, p[0] + 1, p[1] + 1))
        volume = np.zeros(N**2)

        for j in range(N):
            for i in range(N):
                m = i + j * N
                xi[m,...] = (quad_nodes[0][:,np.newaxis].repeat(p[1]+1, axis=1) + 1)\
                          * (space.nodes[0][i+1]-space.nodes[0][i]
                          )/2 + space.nodes[0][i]
                et[m,...] = (quad_nodes[1][np.newaxis,:].repeat(p[0]+1, axis=0) + 1)\
                          * (space.nodes[1][j+1]-space.nodes[1][j]
                          )/2 + space.nodes[1][j]
                volume[m] = (space.nodes[0][i+1]-space.nodes[0][i]) \
                          * (space.nodes[1][j+1]-space.nodes[1][j])  * magic_factor

        self.___D_cache___[N] = xi, et, volume, quad_weights

        return self.___D_cache___[N]



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

        JC = dict()
        LLC = dict() # root-cell-wise Local cochain

        rcfc = mesh.rcfc

        for rp in rcfc:
            cell = rcfc[rp]
            N = cell.N
            xi, et, volume, quad_weights = self.___Pr_prepare_discretize_data___(N)

            xyz = cell.coordinate_transformation.mapping(xi, et)
            typeWr2Metric = cell.___Pr_metric_N_key___

            if typeWr2Metric in JC:
                detJ = JC[typeWr2Metric]
            else:
                detJ = cell.coordinate_transformation.Jacobian(xi, et)
                if not typeWr2Metric.isnumeric():
                    JC[typeWr2Metric] = detJ

            f_xyz = F(*xyz)
            LLC[rp] = np.einsum('jkl, k, l, j -> j',
                f_xyz*detJ, quad_weights[0], quad_weights[1], volume, optimize='greedy')

        return LLC










if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
