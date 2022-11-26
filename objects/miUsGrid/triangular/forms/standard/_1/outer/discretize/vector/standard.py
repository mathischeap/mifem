# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 6:25 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
import numpy as np
from components.quadrature import Quadrature

class miUsTriangular_oS1F_Discretize_StandardVector(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        self._freeze_self_()

    def ___Pr_discretize_preparation___(self, d_='', quad_degree=None):
        """"""
        p = [self._sf_.space.p + 2 for _ in range(2)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
        p_x, p_y = p

        nodes = self._sf_.space.nodes
        self_p = [self._sf_.space.p, self._sf_.space.p]

        edges_size = [nodes[1:] - nodes[:-1] for _ in range(2)]
        cell_nodes = [(0.5 * (edges_size[i][np.newaxis, :]) * (quad_nodes[i][:, np.newaxis] + 1)
            + nodes[:-1]).ravel('F') for i in range(2)]

        if d_ == 'x':
            quad_xi = np.tile(cell_nodes[0], self_p[1] + 1).reshape(
                (p_x + 1, self_p[0] * (self_p[1] + 1)), order='F')
            quad_eta = np.repeat(nodes[np.newaxis, :], self_p[0], axis=0).ravel('F')
            quad_eta = quad_eta[np.newaxis, :].repeat(p_x + 1, axis=0)
            ES = np.tile(edges_size[0], self_p[1] + 1)
            return quad_xi, quad_eta, ES, quad_weights
        elif d_ == 'y':
            quad_xi = np.tile(nodes, self_p[1])[np.newaxis, :].repeat(p_y + 1, axis=0)
            quad_eta = np.repeat(cell_nodes[1].reshape(
                (p_y + 1, self_p[1]), order='F'), self_p[0] + 1, axis=1)
            ES = np.repeat(edges_size[1], self_p[0] + 1)
            return quad_xi, quad_eta, ES, quad_weights
        else:
            raise Exception()

    def __call__(self, update_cochain=True, target='CF', quad_degree=None):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.

        :param update_cochain:
        :param target:
        :param quad_degree:
        :return:
        """
        SELF = self._sf_

        if self.___DISCRETIZE_STANDARD_CACHE___ is None or \
            quad_degree != self.___DISCRETIZE_STANDARD_CACHE___['quadDegree']:
            self.___DISCRETIZE_STANDARD_CACHE___ = dict()

            xi, eta, edge_size_d_xi, quad_weights = \
                self.___Pr_discretize_preparation___(d_='x', quad_degree=quad_degree)
            p = self._sf_.p
            xi = xi[:, p:]                         # we cut the singular dx edges
            eta = eta[:, p:]                       # we cut the singular dx edges
            edge_size_d_xi = edge_size_d_xi[p:]    # we cut the singular dx edges
            self.___DISCRETIZE_STANDARD_CACHE___['X'] = (xi, eta)

            xi, eta, edge_size_d_eta, quad_weights = \
                self.___Pr_discretize_preparation___(d_='y', quad_degree=quad_degree)
            self.___DISCRETIZE_STANDARD_CACHE___['Y'] = (xi, eta)

            edge_size = (edge_size_d_xi, edge_size_d_eta)
            self.___DISCRETIZE_STANDARD_CACHE___['edge'] = edge_size
            self.___DISCRETIZE_STANDARD_CACHE___['quad_weights'] = quad_weights
            self.___DISCRETIZE_STANDARD_CACHE___['quadDegree'] = quad_degree
        else:
            pass

        xi_x, eta_x = self.___DISCRETIZE_STANDARD_CACHE___['X']
        xi_y, eta_y = self.___DISCRETIZE_STANDARD_CACHE___['Y']
        quad_weights = self.___DISCRETIZE_STANDARD_CACHE___['quad_weights']
        edge_size = self.___DISCRETIZE_STANDARD_CACHE___['edge']

        local_dx = dict()
        local_dy = dict()

        # --- target --------------------------------------------------------
        if target == 'CF':
            FUNC = SELF.CF.do.evaluate_func_at_time()
        elif target == 'BC':
            FUNC = SELF.BC.CF.do.evaluate_func_at_time()
        else:
            raise NotImplementedError(f"I cannot deal with target = {target}.")
        # =======================================================================

        for i in SELF.mesh.elements.indices:
            element = SELF.mesh.elements[i]

            smctm = element.coordinate_transformation.mapping(xi_x, eta_x)
            J = element.coordinate_transformation.Jacobian_matrix(xi_x, eta_x)
            J = (J[0][0], J[1][0])
            u = FUNC[0](*smctm)
            v = FUNC[1](*smctm)
            local_dx[i] = np.einsum(
                'jk, j, k -> k', J[0]*v - J[1]*u, quad_weights[0], edge_size[0] * 0.5, optimize='greedy'
            )

            smctm = element.coordinate_transformation.mapping(xi_y, eta_y)
            J = element.coordinate_transformation.Jacobian_matrix(xi_y, eta_y)
            J = (J[0][1], J[1][1])
            u = FUNC[0](*smctm)
            v = FUNC[1](*smctm)
            local_dy[i] = np.einsum(
                'jk, j, k -> k', - J[0]*v + J[1]*u, quad_weights[1], edge_size[1]*0.5, optimize='greedy'
            )

        tM = self._sf_.IDT.transition_matrices

        # give it to cochain.local ...
        cochainLocal = dict()

        for i in SELF.mesh.elements.indices:
            local_cochain = np.hstack((local_dy[i], local_dx[i])) # the difference from inner.
            if i not in tM:
                pass
            else:
                local_cochain = tM[i] @ local_cochain

            cochainLocal[i] = local_cochain

        if update_cochain:
            SELF.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
