

from screws.freeze.base import FrozenOnly
import numpy as np
from screws.quadrature import Quadrature


class _3dCSCG_1Tr_Discretize_StandardVector_T_para(FrozenOnly):
    """"""
    def __init__(self, Tr):
        self._Tr_ = Tr
        self.___cache_DISCRETIZE_STANDARD___ = None
        self._freeze_self_()

    def __call__(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the Trace_parallel component of a standard vector field to all trace
        elements.

        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.
        """

        SELF = self._Tr_

        if target in ('BC',): assert update_cochain is False, \
            f"CANNOT update cochain when target is {target}"

        if self.___cache_DISCRETIZE_STANDARD___ is None or \
            self.___cache_DISCRETIZE_STANDARD___['quadDegree'] != quad_degree:
            p = [SELF.dqp[i] + 1 for i in range(SELF.ndim)] if quad_degree is None else quad_degree
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
            nodes = SELF.space.nodes
            num_edges = [len(nodes[i])-1 for i in range(SELF.ndim)]
            lens = [nodes[i][1:]-nodes[i][0:-1] for i in range(SELF.ndim)]
            qnodes = []
            for i in range(SELF.ndim):
                qnodes_i = ((np.array(quad_nodes[i])+1)/2)[np.newaxis,:].repeat(num_edges[i],
                           axis=0)*lens[i][:,np.newaxis]
                qnodes_i += np.array(nodes[i][:-1])[:,np.newaxis].repeat(p[i]+1, axis=1)
                qnodes.append(qnodes_i)

            # NS ------------------------------------------------------
            qn_NS_dy_y = []
            qn_NS_dy_z = []
            for k in range(SELF.p[2]+1):
                for j in range(SELF.p[1]):
                    qn_NS_dy_y.append(qnodes[1][j])
                    qn_NS_dy_z.append(nodes[2][k]*np.ones(p[1]+1))
            qn_NS_dy_y, qn_NS_dy_z = np.array(qn_NS_dy_y), np.array(qn_NS_dy_z)
            lens_NS_dy = np.tile(lens[1]*0.5, (SELF.p[2]+1))


            qn_NS_dz_y = []
            qn_NS_dz_z = []
            for k in range(SELF.p[2]):
                for j in range(SELF.p[1]+1):
                    qn_NS_dz_y.append(nodes[1][j]*np.ones(p[2]+1))
                    qn_NS_dz_z.append(qnodes[2][k])
            qn_NS_dz_y, qn_NS_dz_z = np.array(qn_NS_dz_y), np.array(qn_NS_dz_z)
            lens_NS_dz = np.repeat(lens[2]*0.5, (SELF.p[1] + 1))

            # WE --------------------------------------------------------------------------------
            qn_WE_dx_x = []
            qn_WE_dx_z = []
            for k in range(SELF.p[2]+1):
                for i in range(SELF.p[0]):
                    qn_WE_dx_x.append(qnodes[0][i])
                    qn_WE_dx_z.append(nodes[2][k]*np.ones(p[0]+1))
            qn_WE_dx_x, qn_WE_dx_z = np.array(qn_WE_dx_x), np.array(qn_WE_dx_z)
            lens_WE_dx = np.tile(lens[0]*0.5, (SELF.p[2] + 1))


            qn_WE_dz_x = []
            qn_WE_dz_z = []
            for k in range(SELF.p[2]):
                for i in range(SELF.p[0]+1):
                    qn_WE_dz_x.append(nodes[0][i]*np.ones(p[2]+1))
                    qn_WE_dz_z.append(qnodes[2][k])
            qn_WE_dz_x, qn_WE_dz_z = np.array(qn_WE_dz_x), np.array(qn_WE_dz_z)
            lens_WE_dz = np.repeat(lens[2]*0.5, (SELF.p[0] + 1))

            #BF --------------------------------------------------------------------------------------------
            qn_BF_dx_x = []
            qn_BF_dx_y = []
            for j in range(SELF.p[1]+1):
                for i in range(SELF.p[0]):
                    qn_BF_dx_x.append(qnodes[0][i])
                    qn_BF_dx_y.append(nodes[1][j]*np.ones(p[0]+1))
            qn_BF_dx_x, qn_BF_dx_y = np.array(qn_BF_dx_x), np.array(qn_BF_dx_y)
            lens_BF_dx = np.tile(lens[0]*0.5, (SELF.p[1] + 1))

            qn_BF_dy_x = []
            qn_BF_dy_y = []
            for j in range(SELF.p[1]):
                for i in range(SELF.p[0]+1):
                    qn_BF_dy_x.append(nodes[0][i]*np.ones(p[1]+1))
                    qn_BF_dy_y.append(qnodes[1][j])
            qn_BF_dy_x, qn_BF_dy_y = np.array(qn_BF_dy_x), np.array(qn_BF_dy_y)
            lens_BF_dy = np.repeat(lens[1]*0.5, (SELF.p[0] + 1))

            LENS = [lens_NS_dy, lens_NS_dz, lens_WE_dx, lens_WE_dz, lens_BF_dx, lens_BF_dy]
            cd = dict()
            cd['quadDegree'] = quad_degree

            cd['LENS'] = LENS
            cd['qn_NS_dy_y'] = qn_NS_dy_y
            cd['qn_NS_dy_z'] = qn_NS_dy_z
            cd['qn_NS_dz_y'] = qn_NS_dz_y
            cd['qn_NS_dz_z'] = qn_NS_dz_z
            cd['qn_WE_dx_x'] = qn_WE_dx_x
            cd['qn_WE_dx_z'] = qn_WE_dx_z
            cd['qn_WE_dz_x'] = qn_WE_dz_x
            cd['qn_WE_dz_z'] = qn_WE_dz_z
            cd['qn_BF_dx_x'] = qn_BF_dx_x
            cd['qn_BF_dx_y'] = qn_BF_dx_y
            cd['qn_BF_dy_x'] = qn_BF_dy_x
            cd['qn_BF_dy_y'] = qn_BF_dy_y
            cd['quad_weights'] = quad_weights

            self.___cache_DISCRETIZE_STANDARD___ = cd

        else:
            LENS = self.___cache_DISCRETIZE_STANDARD___['LENS']
            qn_NS_dy_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_dy_y']
            qn_NS_dy_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_dy_z']
            qn_NS_dz_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_dz_y']
            qn_NS_dz_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_dz_z']
            qn_WE_dx_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_dx_x']
            qn_WE_dx_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_dx_z']
            qn_WE_dz_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_dz_x']
            qn_WE_dz_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_dz_z']
            qn_BF_dx_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_dx_x']
            qn_BF_dx_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_dx_y']
            qn_BF_dy_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_dy_x']
            qn_BF_dy_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_dy_y']
            quad_weights = self.___cache_DISCRETIZE_STANDARD___['quad_weights']

        if target == 'func':
            assert SELF.func.body is not None, f"No func.body!"

        else:
            raise NotImplementedError(f"1Trace = discretize_VectorField_standard: "
                                      f"Not applicable for target={target}.")


        lens_NS_dy, lens_NS_dz, lens_WE_dx, lens_WE_dz, lens_BF_dx, lens_BF_dy = LENS
        local_TEW = dict()
        for key in SELF.mesh.trace.elements:
            te = SELF.mesh.trace.elements[key]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side

            if ele_side in 'NS':

                J = te.coordinate_transformation.Jacobian_matrix(qn_NS_dy_y, qn_NS_dy_z)
                J = (J[0][0], J[1][0], J[2][0]) # dy of (dy, dz)
                x, y, z = te.coordinate_transformation.mapping(qn_NS_dy_y, qn_NS_dy_z,
                                                               from_element=ele, side=ele_side)
                u = SELF.func.body[0](x, y, z)
                v = SELF.func.body[1](x, y, z)
                w = SELF.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[1]
                C = lens_NS_dy
                cochain_dy = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')


                J = te.coordinate_transformation.Jacobian_matrix(qn_NS_dz_y, qn_NS_dz_z)
                J = (J[0][1], J[1][1], J[2][1]) # dz of (dy, dz)
                x, y, z = te.coordinate_transformation.mapping(qn_NS_dz_y, qn_NS_dz_z,
                                                               from_element=ele, side=ele_side)
                u = SELF.func.body[0](x, y, z)
                v = SELF.func.body[1](x, y, z)
                w = SELF.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[2]
                C = lens_NS_dz
                cochain_dz = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dy, cochain_dz])


            elif ele_side in 'WE':

                J = te.coordinate_transformation.Jacobian_matrix(qn_WE_dx_x, qn_WE_dx_z)
                J = (J[0][1], J[1][1], J[2][1]) # dx of (dz, dx)
                x, y, z = te.coordinate_transformation.mapping(qn_WE_dx_x, qn_WE_dx_z,
                                                               from_element=ele, side=ele_side)
                u = SELF.func.body[0](x, y, z)
                v = SELF.func.body[1](x, y, z)
                w = SELF.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[0]
                C = lens_WE_dx
                cochain_dx = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')


                J = te.coordinate_transformation.Jacobian_matrix(qn_WE_dz_x, qn_WE_dz_z)
                J = (J[0][0], J[1][0], J[2][0]) # dz of (dz, dx)
                x, y, z = te.coordinate_transformation.mapping(qn_WE_dz_x, qn_WE_dz_z,
                                                               from_element=ele, side=ele_side)
                u = SELF.func.body[0](x, y, z)
                v = SELF.func.body[1](x, y, z)
                w = SELF.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[2]
                C = lens_WE_dz
                cochain_dz = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dx, cochain_dz])


            elif ele_side in 'BF':

                J = te.coordinate_transformation.Jacobian_matrix(qn_BF_dx_x, qn_BF_dx_y)
                J = (J[0][0], J[1][0], J[2][0]) # dx of (dx, dy)
                x, y, z = te.coordinate_transformation.mapping(qn_BF_dx_x, qn_BF_dx_y,
                                                               from_element=ele, side=ele_side)
                u = SELF.func.body[0](x, y, z)
                v = SELF.func.body[1](x, y, z)
                w = SELF.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[0]
                C = lens_BF_dx
                cochain_dx = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')


                J = te.coordinate_transformation.Jacobian_matrix(qn_BF_dy_x, qn_BF_dy_y)
                J = (J[0][1], J[1][1], J[2][1]) # dy of (dx, dy)
                x, y, z = te.coordinate_transformation.mapping(qn_BF_dy_x, qn_BF_dy_y,
                                                               from_element=ele, side=ele_side)
                u = SELF.func.body[0](x, y, z)
                v = SELF.func.body[1](x, y, z)
                w = SELF.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[1]
                C = lens_BF_dy
                cochain_dy = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dx, cochain_dy])

            else:
                raise Exception()

            if not SELF.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[key] = te_primal_local


        if update_cochain: SELF.cochain.local_TEW = local_TEW
        # 'locally full local TEW cochain': provide cochain.local_TEW and for all dofs on the
        # trace element.
        return 'locally full local TEW cochain', local_TEW