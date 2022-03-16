from screws.freeze.inheriting.frozen_only import FrozenOnly

import numpy as np




class _2dCSCG_S1Fi_Discretize_StandardVector(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func', quad_degree=None):
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
                SELF.___PRIVATE_discretize_preparation___(d_='x', quad_degree=quad_degree)
            self.___DISCRETIZE_STANDARD_CACHE___['X'] = (xi, eta)

            xi, eta, edge_size_d_eta, quad_weights = \
                SELF.___PRIVATE_discretize_preparation___(d_='y', quad_degree=quad_degree)
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
        if target == 'func':
            FUNC = SELF.func.body
        else:
            raise NotImplementedError(f"I cannot deal with target = {target}.")
        # =======================================================================

        JXC, JYC = dict(), dict()
        for i in SELF.mesh.elements.indices:
            element = SELF.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark

            smctm = element.coordinate_transformation.mapping(xi_x, eta_x)
            if typeWr2Metric in JXC:
                J = JXC[typeWr2Metric]
            else:
                J = element.coordinate_transformation.Jacobian_matrix(xi_x, eta_x)
                if isinstance(typeWr2Metric, str):
                    JXC[typeWr2Metric] = J
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                u = FUNC[0](*smctm)
                local_dx[i] = np.einsum(
                    'jk, j, k -> k', J[0][0]*u, quad_weights[0], edge_size[0] * 0.5, optimize='greedy'
                )
            else:
                J = (J[0][0], J[1][0])
                u = FUNC[0](*smctm)
                v = FUNC[1](*smctm)
                local_dx[i] = np.einsum(
                    'jk, j, k -> k', J[0]*u + J[1]*v, quad_weights[0], edge_size[0] * 0.5, optimize='greedy'
                )

            smctm = element.coordinate_transformation.mapping(xi_y, eta_y)
            if typeWr2Metric in JYC:
                J = JYC[typeWr2Metric]
            else:
                J = element.coordinate_transformation.Jacobian_matrix(xi_y, eta_y)
                if isinstance(typeWr2Metric, str):
                    JYC[typeWr2Metric] = J
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                v = FUNC[1](*smctm)
                local_dy[i] = np.einsum(
                    'jk, j, k -> k', J[1][1]*v, quad_weights[1], edge_size[1]*0.5, optimize='greedy'
                )
            else:
                J = (J[0][1], J[1][1])
                u = FUNC[0](*smctm)
                v = FUNC[1](*smctm)
                local_dy[i] = np.einsum(
                    'jk, j, k -> k', J[0]*u + J[1]*v, quad_weights[1], edge_size[1]*0.5, optimize='greedy'
                )

        del JXC, JYC
        # isisKronecker? ...
        if not SELF.space.IS_Kronecker: raise NotImplementedError()
        # give it to cochain.local ...
        cochainLocal = dict()
        for i in SELF.mesh.elements.indices:
            cochainLocal[i] = np.hstack((local_dx[i], local_dy[i]))
        if update_cochain:
            SELF.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal