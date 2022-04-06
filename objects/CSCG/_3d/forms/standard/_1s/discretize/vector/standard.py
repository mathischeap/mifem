from screws.freeze.base import FrozenOnly
import numpy as np

from screws.quadrature import Quadrature

class _3dCSCG_Discretize_Standard(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        self._freeze_self_()

    def __call__(self, update_cochain=True, quad_degree=None):
        """The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.

        :param update_cochain:
        :param quad_degree:
        :return:
        """
        SELF = self._sf_

        if self.___DISCRETIZE_STANDARD_CACHE___ is None or \
            quad_degree != self.___DISCRETIZE_STANDARD_CACHE___['quadDegree']:
            self.___DISCRETIZE_STANDARD_CACHE___ = dict()

            xi, eta, sigma, edge_size_d_xi, quad_weights = \
                self.___PRIVATE_discretize_preparation___(d_='x', quad_degree=quad_degree)
            self.___DISCRETIZE_STANDARD_CACHE___['X'] = (xi, eta, sigma)

            xi, eta, sigma, edge_size_d_eta, quad_weights = \
                self.___PRIVATE_discretize_preparation___(d_='y', quad_degree=quad_degree)
            self.___DISCRETIZE_STANDARD_CACHE___['Y'] = (xi, eta, sigma)

            xi, eta, sigma, edge_size_d_sigma, quad_weights = \
                self.___PRIVATE_discretize_preparation___(d_='z', quad_degree=quad_degree)
            self.___DISCRETIZE_STANDARD_CACHE___['Z'] = (xi, eta, sigma)

            edge_size = (edge_size_d_xi, edge_size_d_eta, edge_size_d_sigma)
            self.___DISCRETIZE_STANDARD_CACHE___['edge'] = edge_size
            self.___DISCRETIZE_STANDARD_CACHE___['quad_weights'] = quad_weights
            self.___DISCRETIZE_STANDARD_CACHE___['quadDegree'] = quad_degree
        else:
            pass

        xi_x, eta_x, sigma_x = self.___DISCRETIZE_STANDARD_CACHE___['X']
        xi_y, eta_y, sigma_y = self.___DISCRETIZE_STANDARD_CACHE___['Y']
        xi_z, eta_z, sigma_z = self.___DISCRETIZE_STANDARD_CACHE___['Z']
        quad_weights = self.___DISCRETIZE_STANDARD_CACHE___['quad_weights']
        edge_size = self.___DISCRETIZE_STANDARD_CACHE___['edge']

        local_dx = dict()
        local_dy = dict()
        local_dz = dict()

        JXC, JYC, JZC = dict(), dict(), dict()
        for i in SELF.mesh.elements.indices:
            element = SELF.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark

            smctm = element.coordinate_transformation.mapping(xi_x, eta_x, sigma_x)
            if typeWr2Metric in JXC:
                J = JXC[typeWr2Metric]
            else:
                J = element.coordinate_transformation.Jacobian_matrix(xi_x, eta_x, sigma_x)
                if isinstance(typeWr2Metric, str):
                    JXC[typeWr2Metric] = J
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                u = SELF.func.body[0](*smctm)
                local_dx[i] = np.einsum(
                    'jk, j, k -> k', J[0][0]*u, quad_weights[0],
                    edge_size[0] * 0.5, optimize='greedy'
                )
            else:
                J = (J[0][0], J[1][0], J[2][0])
                u = SELF.func.body[0](*smctm)
                v = SELF.func.body[1](*smctm)
                w = SELF.func.body[2](*smctm)
                local_dx[i] = np.einsum(
                    'jk, j, k -> k', J[0]*u + J[1]*v + J[2]*w, quad_weights[0],
                    edge_size[0] * 0.5, optimize='greedy'
                )

            smctm = element.coordinate_transformation.mapping(xi_y, eta_y, sigma_y)
            if typeWr2Metric in JYC:
                J = JYC[typeWr2Metric]
            else:
                J = element.coordinate_transformation.Jacobian_matrix(xi_y, eta_y, sigma_y)
                if isinstance(typeWr2Metric, str):
                    JYC[typeWr2Metric] = J
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                v = SELF.func.body[1](*smctm)
                local_dy[i] = np.einsum(
                    'jk, j, k -> k', J[1][1]*v, quad_weights[1],
                    edge_size[1]*0.5, optimize='greedy'
                )
            else:
                J = (J[0][1], J[1][1], J[2][1])
                u = SELF.func.body[0](*smctm)
                v = SELF.func.body[1](*smctm)
                w = SELF.func.body[2](*smctm)
                local_dy[i] = np.einsum(
                    'jk, j, k -> k', J[0]*u + J[1]*v + J[2]*w, quad_weights[1],
                    edge_size[1]*0.5, optimize='greedy'
                )

            smctm = element.coordinate_transformation.mapping(xi_z, eta_z, sigma_z)
            if typeWr2Metric in JZC:
                J = JZC[typeWr2Metric]
            else:
                J = element.coordinate_transformation.Jacobian_matrix(xi_z, eta_z, sigma_z)
                if isinstance(typeWr2Metric, str):
                    JZC[typeWr2Metric] = J
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                w = SELF.func.body[2](*smctm)
                local_dz[i] = np.einsum(
                    'jk, j, k -> k', J[2][2]*w, quad_weights[2],
                    edge_size[2]*0.5, optimize='greedy'
                )
            else:
                J = (J[0][2], J[1][2], J[2][2])
                u = SELF.func.body[0](*smctm)
                v = SELF.func.body[1](*smctm)
                w = SELF.func.body[2](*smctm)
                local_dz[i] = np.einsum(
                    'jk, j, k -> k', J[0]*u + J[1]*v + J[2]*w, quad_weights[2],
                    edge_size[2]*0.5, optimize='greedy'
                )
        del JXC, JYC, JZC
        # isisKronecker? ...
        if not SELF.space.IS_Kronecker: raise NotImplementedError()
        # give it to cochain.local ...
        cochainLocal = dict()
        for i in SELF.mesh.elements.indices:
            cochainLocal[i] = np.hstack((local_dx[i], local_dy[i], local_dz[i]))
        if update_cochain:
            SELF.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal

    def ___PRIVATE_discretize_preparation___(self, d_='', quad_degree=None):
        SELF = self._sf_
        p = [SELF.dqp[i] + 1 for i in range(SELF.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
        quad_num_nodes = [len(quad_nodes_i) for quad_nodes_i in quad_nodes]
        sbn0 = SELF.space.nodes[0]
        sbn1 = SELF.space.nodes[1]
        sbn2 = SELF.space.nodes[2]
        if d_ == 'x':
            a = sbn0[1:] - sbn0[:-1]
            a = a.ravel('F')
            b = (SELF.p[1] + 1) * (SELF.p[2] + 1)
            edge_size_x = np.tile(a, b)
            snbx = b * SELF.p[0]
            D = quad_nodes[0][:, np.newaxis].repeat(snbx, axis=1) + 1
            assert np.shape(D)[1] == len(edge_size_x)
            xi1 = D * edge_size_x / 2
            xi2 = np.tile(sbn0[:-1], b)
            xi = xi1 + xi2
            eta = np.tile(np.tile(sbn1[:, np.newaxis].repeat(quad_num_nodes[0], axis=1).T,
                                  (SELF.p[0], 1)).reshape((quad_num_nodes[0], SELF.p[0] * (SELF.p[1] + 1)),
                                                          order='F'),
                          (1, SELF.p[2] + 1))
            sigma = sbn2.repeat(SELF.p[0] * (SELF.p[1] + 1))[np.newaxis, :].repeat(
                quad_num_nodes[0], axis=0)
            return xi, eta, sigma, edge_size_x, quad_weights
        elif d_ == 'y':
            edge_size_y = np.tile(np.repeat((sbn1[1:] - sbn1[:-1]),
                                            SELF.p[0] + 1), SELF.p[2] + 1)
            xi = np.tile(sbn0, SELF.p[1] * (SELF.p[2] + 1))[np.newaxis, :].repeat(
                quad_num_nodes[1], axis=0)
            sn_by = SELF.num.basis_components[1]
            eta1 = (quad_nodes[1][:, np.newaxis].repeat(sn_by, axis=1) + 1) * edge_size_y / 2
            eta2 = np.tile(np.repeat(sbn1[:-1], (SELF.p[0] + 1)), (SELF.p[2] + 1))
            eta = eta1 + eta2
            sigma = sbn2.repeat(SELF.p[1] * (SELF.p[0] + 1))[np.newaxis, :].repeat(
                quad_num_nodes[1], axis=0)
            return xi, eta, sigma, edge_size_y, quad_weights
        elif d_ == 'z':
            edge_size_z = np.repeat((sbn2[1:] - sbn2[:-1]),
                                    SELF.p[0] + 1).repeat(SELF.p[1] + 1)
            xi = np.tile(sbn0, (SELF.p[1] + 1) * (SELF.p[2]))[np.newaxis, :].repeat(
                quad_num_nodes[2], axis=0)
            eta = np.tile(np.repeat(sbn1, (SELF.p[0] + 1)), SELF.p[2])[np.newaxis, :].repeat(
                quad_num_nodes[2], axis=0)
            sn_bz = SELF.num.basis_components[2]
            sigma1 = (quad_nodes[2][:, np.newaxis].repeat(sn_bz, axis=1) + 1) * edge_size_z / 2
            sigma2 = sbn2[:-1].repeat((SELF.p[0] + 1) * (SELF.p[1] + 1))
            sigma = sigma1 + sigma2
            return xi, eta, sigma, edge_size_z, quad_weights
        else:
            raise Exception()