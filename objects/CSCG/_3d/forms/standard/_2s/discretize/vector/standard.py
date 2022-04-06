from screws.freeze.base import FrozenOnly
import numpy as np

from screws.quadrature import Quadrature



class _3dCSCG_Discretize_Standard(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        self._freeze_self_()


    def __call__(self, update_cochain=True, target='func', quad_degree=None):
        """The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.

        :param update_cochain:
        :param target:
        :param quad_degree:
        :return:
        """
        SELF = self._sf_
        p = [SELF.dqp[i] + 1 for i in range(SELF.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None or \
                quad_degree != self.___DISCRETIZE_STANDARD_CACHE___['quadDegree']:
            self.___DISCRETIZE_STANDARD_CACHE___ = dict()

            xi = np.zeros((SELF.num.basis_components[0], p[1] + 1, p[2] + 1))
            et = np.zeros((SELF.num.basis_components[0], p[1] + 1, p[2] + 1))
            si = np.zeros((SELF.num.basis_components[0], p[1] + 1, p[2] + 1))
            area_dydz = np.zeros((SELF.num.basis_components[0]))
            for k in range(SELF.p[2]):
                for j in range(SELF.p[1]):
                    for i in range(SELF.p[0] + 1):
                        m = i + j * (SELF.p[0] + 1) + k * (SELF.p[0] + 1) * SELF.p[1]
                        xi[m, ...] = np.ones((p[1] + 1, p[2] + 1)) * SELF.space.nodes[0][i]
                        et[m, ...] = (quad_nodes[1][:, np.newaxis].repeat(p[2] + 1, axis=1) + 1) * (
                                SELF.space.nodes[1][j + 1] - SELF.space.nodes[1][j]) / 2 + \
                                     SELF.space.nodes[1][j]
                        si[m, ...] = (quad_nodes[2][np.newaxis, :].repeat((p[1] + 1), axis=0) + 1) * (
                                SELF.space.nodes[2][k + 1] - SELF.space.nodes[2][k]) / 2 + \
                                     SELF.space.nodes[2][k]
                        area_dydz[m] = (SELF.space.nodes[2][k + 1] - SELF.space.nodes[2][k]) \
                                       * (SELF.space.nodes[1][j + 1] - SELF.space.nodes[1][j])

            self.___DISCRETIZE_STANDARD_CACHE___['X'] = (xi, et, si, area_dydz)
            #
            xi = np.zeros((SELF.num.basis_components[1], p[0] + 1, p[2] + 1))
            et = np.zeros((SELF.num.basis_components[1], p[0] + 1, p[2] + 1))
            si = np.zeros((SELF.num.basis_components[1], p[0] + 1, p[2] + 1))
            area_dzdx = np.zeros((SELF.num.basis_components[1]))
            for k in range(SELF.p[2]):
                for j in range(SELF.p[1] + 1):
                    for i in range(SELF.p[0]):
                        m = i + j * SELF.p[0] + k * (SELF.p[1] + 1) * SELF.p[0]
                        xi[m, ...] = (quad_nodes[0][:, np.newaxis].repeat(p[2] + 1, axis=1) + 1) * (
                                SELF.space.nodes[0][i + 1] - SELF.space.nodes[0][i]) / 2 + \
                                     SELF.space.nodes[0][i]
                        et[m, ...] = np.ones((p[0] + 1, p[2] + 1)) * SELF.space.nodes[1][j]
                        si[m, ...] = (quad_nodes[2][np.newaxis, :].repeat(p[0] + 1, axis=0) + 1) * (
                                SELF.space.nodes[2][k + 1] - SELF.space.nodes[2][k]) / 2 + \
                                     SELF.space.nodes[2][k]
                        area_dzdx[m] = (SELF.space.nodes[2][k + 1] - SELF.space.nodes[2][k]) \
                                       * (SELF.space.nodes[0][i + 1] - SELF.space.nodes[0][i])

            self.___DISCRETIZE_STANDARD_CACHE___['Y'] = (xi, et, si, area_dzdx)

            xi = np.zeros((SELF.num.basis_components[2], p[0] + 1, p[1] + 1))
            et = np.zeros((SELF.num.basis_components[2], p[0] + 1, p[1] + 1))
            si = np.zeros((SELF.num.basis_components[2], p[0] + 1, p[1] + 1))
            area_dxdy = np.zeros((SELF.num.basis_components[2]))
            for k in range(SELF.p[2] + 1):
                for j in range(SELF.p[1]):
                    for i in range(SELF.p[0]):
                        m = i + j * SELF.p[0] + k * SELF.p[1] * SELF.p[0]
                        xi[m, ...] = (quad_nodes[0][:, np.newaxis].repeat(p[1] + 1, axis=1) + 1) * (
                                SELF.space.nodes[0][i + 1] - SELF.space.nodes[0][i]) / 2 + \
                                     SELF.space.nodes[0][i]
                        et[m, ...] = (quad_nodes[1][np.newaxis, :].repeat(p[0] + 1, axis=0) + 1) * (
                                SELF.space.nodes[1][j + 1] - SELF.space.nodes[1][j]) / 2 + \
                                     SELF.space.nodes[1][j]
                        si[m, ...] = np.ones((p[0] + 1, p[1] + 1)) * SELF.space.nodes[2][k]
                        area_dxdy[m] = (SELF.space.nodes[1][j + 1] - SELF.space.nodes[1][j]) \
                                       * (SELF.space.nodes[0][i + 1] - SELF.space.nodes[0][i])

            self.___DISCRETIZE_STANDARD_CACHE___['Z'] = (xi, et, si, area_dxdy)
            self.___DISCRETIZE_STANDARD_CACHE___['quadDegree'] = quad_degree
            del xi, et, si
        else:
            pass

        xi_x, et_x, si_x, area_dydz = self.___DISCRETIZE_STANDARD_CACHE___['X']
        xi_y, et_y, si_y, area_dzdx = self.___DISCRETIZE_STANDARD_CACHE___['Y']
        xi_z, et_z, si_z, area_dxdy = self.___DISCRETIZE_STANDARD_CACHE___['Z']

        local_dydz = dict()
        local_dzdx = dict()
        local_dxdy = dict()

        JXC, JYC, JZC = dict(), dict(), dict()

        if target == 'func':
            FUNC = SELF.func
        elif target == 'BC':
            FUNC = SELF.BC
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"_2Form.___PRIVATE_discretize_standard_ftype___ "
                f"does not work for target={target}.")

        for i in SELF.mesh.elements:
            element = SELF.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark

            smctm_x = element.coordinate_transformation.mapping(xi_x, et_x, si_x)
            if typeWr2Metric in JXC:
                Jx_0, Jx_1, Jx_2 = JXC[typeWr2Metric]
            else:
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    J11 = element.coordinate_transformation.J11(xi_x, et_x, si_x)
                    J22 = element.coordinate_transformation.J22(xi_x, et_x, si_x)
                    Jx_0 = J11 * J22
                    Jx_1 = 0
                    Jx_2 = 0
                else:
                    J = element.coordinate_transformation.Jacobian_matrix(xi_x, et_x, si_x)
                    Jx_0 = J[1][1] * J[2][2] - J[1][2] * J[2][1]
                    Jx_1 = J[2][1] * J[0][2] - J[2][2] * J[0][1]
                    Jx_2 = J[0][1] * J[1][2] - J[0][2] * J[1][1]
                if isinstance(typeWr2Metric, str):
                    JXC[typeWr2Metric] = Jx_0, Jx_1, Jx_2

            smctm_y = element.coordinate_transformation.mapping(xi_y, et_y, si_y)
            if typeWr2Metric in JYC:
                Jy_0, Jy_1, Jy_2 = JYC[typeWr2Metric]
            else:
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    J00 = element.coordinate_transformation.J00(xi_y, et_y, si_y)
                    J22 = element.coordinate_transformation.J22(xi_y, et_y, si_y)
                    Jy_0 = 0
                    Jy_1 = J00 * J22
                    Jy_2 = 0
                else:
                    J = element.coordinate_transformation.Jacobian_matrix(xi_y, et_y, si_y)
                    Jy_0 = J[1][2] * J[2][0] - J[1][0] * J[2][2]
                    Jy_1 = J[2][2] * J[0][0] - J[2][0] * J[0][2]
                    Jy_2 = J[0][2] * J[1][0] - J[0][0] * J[1][2]
                if isinstance(typeWr2Metric, str):
                    JYC[typeWr2Metric] = Jy_0, Jy_1, Jy_2

            smctm_z = element.coordinate_transformation.mapping(xi_z, et_z, si_z)
            if typeWr2Metric in JZC:
                Jz_0, Jz_1, Jz_2 = JZC[typeWr2Metric]
            else:
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    J11 = element.coordinate_transformation.J11(xi_z, et_z, si_z)
                    J00 = element.coordinate_transformation.J00(xi_z, et_z, si_z)
                    Jz_0 = 0
                    Jz_1 = 0
                    Jz_2 = J11 * J00
                else:
                    J = element.coordinate_transformation.Jacobian_matrix(xi_z, et_z, si_z)
                    Jz_0 = J[1][0] * J[2][1] - J[1][1] * J[2][0]
                    Jz_1 = J[2][0] * J[0][1] - J[2][1] * J[0][0]
                    Jz_2 = J[0][0] * J[1][1] - J[0][1] * J[1][0]
                if isinstance(typeWr2Metric, str):
                    JZC[typeWr2Metric] = Jz_0, Jz_1, Jz_2

            # Now it is time to do the reduction.

            u = FUNC.body[0](*smctm_x)
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                uvw_dydz = Jx_0 * u
            else:
                v = FUNC.body[1](*smctm_x)
                w = FUNC.body[2](*smctm_x)
                uvw_dydz = Jx_0 * u + Jx_1 * v + Jx_2 * w
            local_dydz[i] = self.___PRIVATE_discretize_standard_einsum___(
                uvw_dydz, quad_weights[1], quad_weights[2], area_dydz
            )
            v = FUNC.body[1](*smctm_y)
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                uvw_dzdx = Jy_1 * v
            else:
                u = FUNC.body[0](*smctm_y)
                w = FUNC.body[2](*smctm_y)
                uvw_dzdx = Jy_0 * u + Jy_1 * v + Jy_2 * w
            local_dzdx[i] = self.___PRIVATE_discretize_standard_einsum___(
                uvw_dzdx, quad_weights[0], quad_weights[2], area_dzdx
            )
            w = FUNC.body[2](*smctm_z)
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                uvw_dxdy = Jz_2 * w
            else:
                u = FUNC.body[0](*smctm_z)
                v = FUNC.body[1](*smctm_z)
                uvw_dxdy = Jz_0 * u + Jz_1 * v + Jz_2 * w
            local_dxdy[i] = self.___PRIVATE_discretize_standard_einsum___(
                uvw_dxdy, quad_weights[0], quad_weights[1], area_dxdy
            )
        del JXC, JYC, JZC
        # isKronecker? ...
        if not SELF.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        cochainLocal = dict()
        for i in SELF.mesh.elements.indices:
            cochainLocal[i] = np.hstack((local_dydz[i], local_dzdx[i], local_dxdy[i]))
        if update_cochain: SELF.cochain.local = cochainLocal
        # 'locally full local cochain': provide cochain.local and locally they are full for all local dofs
        return 'locally full local cochain', cochainLocal


    @staticmethod
    def ___PRIVATE_discretize_standard_einsum___(uvw, quad_weights_1, quad_weights_2, area):
        """ """
        return np.einsum('jkl, kl, j -> j',
                         uvw, np.tensordot(quad_weights_1, quad_weights_2, axes=0),
                         area * 0.25, optimize='optimal'
                         )