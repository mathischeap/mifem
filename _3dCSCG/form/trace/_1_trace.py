# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
from abc import ABC
if './' not in sys.path: sys.path.append('./')

from root.config import *
from SCREWS.quadrature import Quadrature
from _3dCSCG.form.trace.main import _3dCSCG_Standard_Trace
from scipy import sparse as spspa

class _1Trace(_3dCSCG_Standard_Trace, ABC):
    """
    Trace 1-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-1-trace-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_trace_1form')
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self.___cache_DISCRETIZE_STANDARD___ = None
        self.___cache_DISCRETIZE_TEW___ = None
        super().RESET_cache()

    def ___TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard', 'trace-element-wise'), \
                f"3dCSCG 1-trace FUNC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."

        else:
            raise NotImplementedError(
                f"1-trace form cannot accommodate {func_body}.")

    def discretize(self, update_cochain=True, target='func', component='T_para', **kwargs):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param target:
        :param component:
        :param kwargs: Keywords arguments to be passed to particular discretization schemes.
        :return: The cochain corresponding to the particular discretization scheme.
        """



        if target == 'func':

            if self.TW.func.body.__class__.__name__ == '_3dCSCG_VectorField':


                if self.func.ftype == 'standard':

                    if component == 'T_para':

                        return self.___PRIVATE_discretize_VectorField_of_ftype_standard_at_component_T_para___(
                            update_cochain=update_cochain, target='func', **kwargs)

                    elif component == 'T_perp':

                        return self.___PRIVATE_discretize_VectorField_of_ftype_standard_at_component_T_perp___(
                            update_cochain=update_cochain, target='func', **kwargs)

                    else:
                        raise Exception(f"1-trace-form discretization targeting vector func of standard type. "
                                        f"I cannot discretize component = {component}. It should be either"
                                        f"'T_para' (parallel trace) or 'T_perp' (perpendicular trace).")


                elif self.func.ftype == 'trace-element-wise':
                    # we do not care this trace-element-wise vector is T_para or T_perp vector, we just discretize it to the trace.
                    return self.___PRIVATE_discretize_VectorField_of_ftype_trace_element_wise___(
                        update_cochain=update_cochain, target='func', **kwargs)


                else:
                    raise Exception(f'3dCSCG 1-trace can not (target func) discretize _3dCSCG_VectorField of ftype {self.func.ftype}.')



            else:
                raise NotImplementedError(f'3dCSCG 1-trace can not (target func) discretize {self.TW.func.body.__class__}.')

        else:
            raise NotImplementedError(f"target={target} not implemented "
                                      f"for 3d CSCG 1-trace form discretization.")


    def ___PRIVATE_discretize_VectorField_of_ftype_standard_at_component_T_para___(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the a the Trace_parallel component of a standard vector field to all trace
        elements.

        """
        if target in ('BC',): assert update_cochain is False, f"CANNOT update cochain when target is {target}"

        if self.___cache_DISCRETIZE_STANDARD___ is None or \
            self.___cache_DISCRETIZE_STANDARD___['quadDegree'] != quad_degree:
            p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
            nodes = self.space.nodes
            num_edges = [len(nodes[i])-1 for i in range(self.ndim)]
            lens = [nodes[i][1:]-nodes[i][0:-1] for i in range(self.ndim)]
            qnodes = []
            for i in range(self.ndim):
                qnodes_i = ((np.array(quad_nodes[i])+1)/2)[np.newaxis,:].repeat(num_edges[i],
                           axis=0)*lens[i][:,np.newaxis]
                qnodes_i += np.array(nodes[i][:-1])[:,np.newaxis].repeat(p[i]+1, axis=1)
                qnodes.append(qnodes_i)

            # NS ------------------------------------------------------
            qn_NS_dy_y = []
            qn_NS_dy_z = []
            for k in range(self.p[2]+1):
                for j in range(self.p[1]):
                    qn_NS_dy_y.append(qnodes[1][j])
                    qn_NS_dy_z.append(nodes[2][k]*np.ones(p[1]+1))
            qn_NS_dy_y, qn_NS_dy_z = np.array(qn_NS_dy_y), np.array(qn_NS_dy_z)
            lens_NS_dy = np.tile(lens[1]*0.5, (self.p[2]+1))


            qn_NS_dz_y = []
            qn_NS_dz_z = []
            for k in range(self.p[2]):
                for j in range(self.p[1]+1):
                    qn_NS_dz_y.append(nodes[1][j]*np.ones(p[2]+1))
                    qn_NS_dz_z.append(qnodes[2][k])
            qn_NS_dz_y, qn_NS_dz_z = np.array(qn_NS_dz_y), np.array(qn_NS_dz_z)
            lens_NS_dz = np.repeat(lens[2]*0.5, (self.p[1] + 1))

            # WE --------------------------------------------------------------------------------
            qn_WE_dx_x = []
            qn_WE_dx_z = []
            for k in range(self.p[2]+1):
                for i in range(self.p[0]):
                    qn_WE_dx_x.append(qnodes[0][i])
                    qn_WE_dx_z.append(nodes[2][k]*np.ones(p[0]+1))
            qn_WE_dx_x, qn_WE_dx_z = np.array(qn_WE_dx_x), np.array(qn_WE_dx_z)
            lens_WE_dx = np.tile(lens[0]*0.5, (self.p[2] + 1))


            qn_WE_dz_x = []
            qn_WE_dz_z = []
            for k in range(self.p[2]):
                for i in range(self.p[0]+1):
                    qn_WE_dz_x.append(nodes[0][i]*np.ones(p[2]+1))
                    qn_WE_dz_z.append(qnodes[2][k])
            qn_WE_dz_x, qn_WE_dz_z = np.array(qn_WE_dz_x), np.array(qn_WE_dz_z)
            lens_WE_dz = np.repeat(lens[2]*0.5, (self.p[0] + 1))

            #BF --------------------------------------------------------------------------------------------
            qn_BF_dx_x = []
            qn_BF_dx_y = []
            for j in range(self.p[1]+1):
                for i in range(self.p[0]):
                    qn_BF_dx_x.append(qnodes[0][i])
                    qn_BF_dx_y.append(nodes[1][j]*np.ones(p[0]+1))
            qn_BF_dx_x, qn_BF_dx_y = np.array(qn_BF_dx_x), np.array(qn_BF_dx_y)
            lens_BF_dx = np.tile(lens[0]*0.5, (self.p[1] + 1))

            qn_BF_dy_x = []
            qn_BF_dy_y = []
            for j in range(self.p[1]):
                for i in range(self.p[0]+1):
                    qn_BF_dy_x.append(nodes[0][i]*np.ones(p[1]+1))
                    qn_BF_dy_y.append(qnodes[1][j])
            qn_BF_dy_x, qn_BF_dy_y = np.array(qn_BF_dy_x), np.array(qn_BF_dy_y)
            lens_BF_dy = np.repeat(lens[1]*0.5, (self.p[0] + 1))

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
            assert self.func.body is not None, f"No func.body!"
        else:
            raise NotImplementedError(f"1Trace = discretize_VectorField_standard: "
                                      f"Not applicable for target={target}.")

        lens_NS_dy, lens_NS_dz, lens_WE_dx, lens_WE_dz, lens_BF_dx, lens_BF_dy = LENS
        local_TEW = dict()
        for key in self.mesh.trace.elements:
            te = self.mesh.trace.elements[key]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side

            if ele_side in 'NS':

                J = te.coordinate_transformation.Jacobian_matrix(qn_NS_dy_y, qn_NS_dy_z)
                J = (J[0][0], J[1][0], J[2][0]) # dy of (dy, dz)
                x, y, z = te.coordinate_transformation.mapping(qn_NS_dy_y, qn_NS_dy_z, from_element=ele, side=ele_side)
                u = self.func.body[0](x, y, z)
                v = self.func.body[1](x, y, z)
                w = self.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[1]
                C = lens_NS_dy
                cochain_dy = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')


                J = te.coordinate_transformation.Jacobian_matrix(qn_NS_dz_y, qn_NS_dz_z)
                J = (J[0][1], J[1][1], J[2][1]) # dz of (dy, dz)
                x, y, z = te.coordinate_transformation.mapping(qn_NS_dz_y, qn_NS_dz_z, from_element=ele, side=ele_side)
                u = self.func.body[0](x, y, z)
                v = self.func.body[1](x, y, z)
                w = self.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[2]
                C = lens_NS_dz
                cochain_dz = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dy, cochain_dz])


            elif ele_side in 'WE':

                J = te.coordinate_transformation.Jacobian_matrix(qn_WE_dx_x, qn_WE_dx_z)
                J = (J[0][1], J[1][1], J[2][1]) # dx of (dz, dx)
                x, y, z = te.coordinate_transformation.mapping(qn_WE_dx_x, qn_WE_dx_z, from_element=ele, side=ele_side)
                u = self.func.body[0](x, y, z)
                v = self.func.body[1](x, y, z)
                w = self.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[0]
                C = lens_WE_dx
                cochain_dx = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')


                J = te.coordinate_transformation.Jacobian_matrix(qn_WE_dz_x, qn_WE_dz_z)
                J = (J[0][0], J[1][0], J[2][0]) # dz of (dz, dx)
                x, y, z = te.coordinate_transformation.mapping(qn_WE_dz_x, qn_WE_dz_z, from_element=ele, side=ele_side)
                u = self.func.body[0](x, y, z)
                v = self.func.body[1](x, y, z)
                w = self.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[2]
                C = lens_WE_dz
                cochain_dz = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dx, cochain_dz])


            elif ele_side in 'BF':

                J = te.coordinate_transformation.Jacobian_matrix(qn_BF_dx_x, qn_BF_dx_y)
                J = (J[0][0], J[1][0], J[2][0]) # dx of (dx, dy)
                x, y, z = te.coordinate_transformation.mapping(qn_BF_dx_x, qn_BF_dx_y, from_element=ele, side=ele_side)
                u = self.func.body[0](x, y, z)
                v = self.func.body[1](x, y, z)
                w = self.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[0]
                C = lens_BF_dx
                cochain_dx = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')


                J = te.coordinate_transformation.Jacobian_matrix(qn_BF_dy_x, qn_BF_dy_y)
                J = (J[0][1], J[1][1], J[2][1]) # dy of (dx, dy)
                x, y, z = te.coordinate_transformation.mapping(qn_BF_dy_x, qn_BF_dy_y, from_element=ele, side=ele_side)
                u = self.func.body[0](x, y, z)
                v = self.func.body[1](x, y, z)
                w = self.func.body[2](x, y, z)
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[1]
                C = lens_BF_dy
                cochain_dy = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dx, cochain_dy])

            else:
                raise Exception()

            if not self.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[key] = te_primal_local


        if update_cochain: self.cochain.local_TEW = local_TEW
        # 'locally full local TEW cochain': provide cochain.local_TEW and for all dofs on the trace element.
        return 'locally full local TEW cochain', local_TEW

    def ___PRIVATE_discretize_VectorField_of_ftype_standard_at_component_T_perp___(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the a the Trace_perpendicular component of a standard vector field to all trace
        elements.

        """
        if target in ('BC',): assert update_cochain is False, f"CANNOT update cochain when target is {target}"
        raise NotImplementedError()

    def ___PRIVATE_discretize_VectorField_of_ftype_trace_element_wise___(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the a the Trace_parallel component of a standard vector field to all trace
        elements.

        """
        # first check `target` and `update_cochain` inputs----------------------------------------
        if target in ('BC',):
            assert update_cochain is False, f"CANNOT update cochain when target is {target}"

        # ----- 3D: prepare and cache or read from cache the data for the numerical integration ----------------
        if self.___cache_DISCRETIZE_STANDARD___ is None or self.___cache_DISCRETIZE_STANDARD___['quadDegree'] != quad_degree:
            p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
            nodes = self.space.nodes
            num_edges = [len(nodes[i])-1 for i in range(self.ndim)]
            lens = [nodes[i][1:]-nodes[i][0:-1] for i in range(self.ndim)]
            qnodes = []
            for i in range(self.ndim):
                qnodes_i = ((np.array(quad_nodes[i])+1)/2)[np.newaxis,:].repeat(num_edges[i],
                           axis=0)*lens[i][:,np.newaxis]
                qnodes_i += np.array(nodes[i][:-1])[:,np.newaxis].repeat(p[i]+1, axis=1)
                qnodes.append(qnodes_i)

            # NS ------------------------------------------------------
            qn_NS_dy_y = []
            qn_NS_dy_z = []
            for k in range(self.p[2]+1):
                for j in range(self.p[1]):
                    qn_NS_dy_y.append(qnodes[1][j])
                    qn_NS_dy_z.append(nodes[2][k]*np.ones(p[1]+1))
            qn_NS_dy_y, qn_NS_dy_z = np.array(qn_NS_dy_y), np.array(qn_NS_dy_z)
            lens_NS_dy = np.tile(lens[1]*0.5, (self.p[2]+1))

            qn_NS_dz_y = []
            qn_NS_dz_z = []
            for k in range(self.p[2]):
                for j in range(self.p[1]+1):
                    qn_NS_dz_y.append(nodes[1][j]*np.ones(p[2]+1))
                    qn_NS_dz_z.append(qnodes[2][k])
            qn_NS_dz_y, qn_NS_dz_z = np.array(qn_NS_dz_y), np.array(qn_NS_dz_z)
            lens_NS_dz = np.repeat(lens[2]*0.5, (self.p[1] + 1))

            # WE --------------------------------------------------------------------------------
            qn_WE_dx_x = []
            qn_WE_dx_z = []
            for k in range(self.p[2]+1):
                for i in range(self.p[0]):
                    qn_WE_dx_x.append(qnodes[0][i])
                    qn_WE_dx_z.append(nodes[2][k]*np.ones(p[0]+1))
            qn_WE_dx_x, qn_WE_dx_z = np.array(qn_WE_dx_x), np.array(qn_WE_dx_z)
            lens_WE_dx = np.tile(lens[0]*0.5, (self.p[2] + 1))


            qn_WE_dz_x = []
            qn_WE_dz_z = []
            for k in range(self.p[2]):
                for i in range(self.p[0]+1):
                    qn_WE_dz_x.append(nodes[0][i]*np.ones(p[2]+1))
                    qn_WE_dz_z.append(qnodes[2][k])
            qn_WE_dz_x, qn_WE_dz_z = np.array(qn_WE_dz_x), np.array(qn_WE_dz_z)
            lens_WE_dz = np.repeat(lens[2]*0.5, (self.p[0] + 1))

            #BF --------------------------------------------------------------------------------------------
            qn_BF_dx_x = []
            qn_BF_dx_y = []
            for j in range(self.p[1]+1):
                for i in range(self.p[0]):
                    qn_BF_dx_x.append(qnodes[0][i])
                    qn_BF_dx_y.append(nodes[1][j]*np.ones(p[0]+1))
            qn_BF_dx_x, qn_BF_dx_y = np.array(qn_BF_dx_x), np.array(qn_BF_dx_y)
            lens_BF_dx = np.tile(lens[0]*0.5, (self.p[1] + 1))

            qn_BF_dy_x = []
            qn_BF_dy_y = []
            for j in range(self.p[1]):
                for i in range(self.p[0]+1):
                    qn_BF_dy_x.append(nodes[0][i]*np.ones(p[1]+1))
                    qn_BF_dy_y.append(qnodes[1][j])
            qn_BF_dy_x, qn_BF_dy_y = np.array(qn_BF_dy_x), np.array(qn_BF_dy_y)
            lens_BF_dy = np.repeat(lens[1]*0.5, (self.p[0] + 1))

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

        # ----- 1D: prepare and cache or read from cache the data for the numerical integration ----------------
        if self.___cache_DISCRETIZE_TEW___ is None or self.___cache_DISCRETIZE_TEW___['quadDegree'] != quad_degree:
            p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
            nodes = self.space.nodes
            num_edges = [len(nodes[i])-1 for i in range(self.ndim)]
            lens = [nodes[i][1:]-nodes[i][0:-1] for i in range(self.ndim)]
            qnodes = []
            for i in range(self.ndim):
                qnodes_i = ((np.array(quad_nodes[i])+1)/2)[np.newaxis,:].repeat(num_edges[i],
                           axis=0)*lens[i][:,np.newaxis]
                qnodes_i += np.array(nodes[i][:-1])[:,np.newaxis].repeat(p[i]+1, axis=1)
                qnodes.append(list(qnodes_i))

            # NS ------------------------------------------------------------------------------------
            qn_NS_dy_y1 = []
            qn_NS_dy_z1 = []
            for k in range(self.p[2]+1):
                for j in range(self.p[1]):
                    qn_NS_dy_y1.append(qnodes[1][j])
                    qn_NS_dy_z1.append([nodes[2][k],])
            lens_NS_dy = np.tile(lens[1]*0.5, (self.p[2]+1))

            qn_NS_dz_y1 = []
            qn_NS_dz_z1 = []
            for k in range(self.p[2]):
                for j in range(self.p[1]+1):
                    qn_NS_dz_y1.append([nodes[1][j],])
                    qn_NS_dz_z1.append(qnodes[2][k])
            lens_NS_dz = np.repeat(lens[2]*0.5, (self.p[1] + 1))

            # WE --------------------------------------------------------------------------------
            qn_WE_dx_x1 = []
            qn_WE_dx_z1 = []
            for k in range(self.p[2]+1):
                for i in range(self.p[0]):
                    qn_WE_dx_x1.append(qnodes[0][i])
                    qn_WE_dx_z1.append([nodes[2][k],])
            lens_WE_dx = np.tile(lens[0]*0.5, (self.p[2] + 1))


            qn_WE_dz_x1 = []
            qn_WE_dz_z1 = []
            for k in range(self.p[2]):
                for i in range(self.p[0]+1):
                    qn_WE_dz_x1.append([nodes[0][i],])
                    qn_WE_dz_z1.append(qnodes[2][k])
            lens_WE_dz = np.repeat(lens[2]*0.5, (self.p[0] + 1))

            #BF --------------------------------------------------------------------------------------------
            qn_BF_dx_x1 = []
            qn_BF_dx_y1 = []
            for j in range(self.p[1]+1):
                for i in range(self.p[0]):
                    qn_BF_dx_x1.append(qnodes[0][i])
                    qn_BF_dx_y1.append([nodes[1][j],])
            lens_BF_dx = np.tile(lens[0]*0.5, (self.p[1] + 1))

            qn_BF_dy_x1 = []
            qn_BF_dy_y1 = []
            for j in range(self.p[1]):
                for i in range(self.p[0]+1):
                    qn_BF_dy_x1.append([nodes[0][i],])
                    qn_BF_dy_y1.append(qnodes[1][j])
            lens_BF_dy = np.repeat(lens[1]*0.5, (self.p[0] + 1))

            LENS = [lens_NS_dy, lens_NS_dz, lens_WE_dx, lens_WE_dz, lens_BF_dx, lens_BF_dy]
            cd = dict()
            cd['quadDegree'] = quad_degree

            cd['LENS'] = LENS
            cd['qn_NS_dy_y1'] = qn_NS_dy_y1
            cd['qn_NS_dy_z1'] = qn_NS_dy_z1
            cd['qn_NS_dz_y1'] = qn_NS_dz_y1
            cd['qn_NS_dz_z1'] = qn_NS_dz_z1
            cd['qn_WE_dx_x1'] = qn_WE_dx_x1
            cd['qn_WE_dx_z1'] = qn_WE_dx_z1
            cd['qn_WE_dz_x1'] = qn_WE_dz_x1
            cd['qn_WE_dz_z1'] = qn_WE_dz_z1
            cd['qn_BF_dx_x1'] = qn_BF_dx_x1
            cd['qn_BF_dx_y1'] = qn_BF_dx_y1
            cd['qn_BF_dy_x1'] = qn_BF_dy_x1
            cd['qn_BF_dy_y1'] = qn_BF_dy_y1
            cd['quad_weights'] = quad_weights
            self.___cache_DISCRETIZE_TEW___ = cd
        else:
            LENS = self.___cache_DISCRETIZE_TEW___['LENS']
            qn_NS_dy_y1 = self.___cache_DISCRETIZE_TEW___['qn_NS_dy_y1']
            qn_NS_dy_z1 = self.___cache_DISCRETIZE_TEW___['qn_NS_dy_z1']
            qn_NS_dz_y1 = self.___cache_DISCRETIZE_TEW___['qn_NS_dz_y1']
            qn_NS_dz_z1 = self.___cache_DISCRETIZE_TEW___['qn_NS_dz_z1']
            qn_WE_dx_x1 = self.___cache_DISCRETIZE_TEW___['qn_WE_dx_x1']
            qn_WE_dx_z1 = self.___cache_DISCRETIZE_TEW___['qn_WE_dx_z1']
            qn_WE_dz_x1 = self.___cache_DISCRETIZE_TEW___['qn_WE_dz_x1']
            qn_WE_dz_z1 = self.___cache_DISCRETIZE_TEW___['qn_WE_dz_z1']
            qn_BF_dx_x1 = self.___cache_DISCRETIZE_TEW___['qn_BF_dx_x1']
            qn_BF_dx_y1 = self.___cache_DISCRETIZE_TEW___['qn_BF_dx_y1']
            qn_BF_dy_x1 = self.___cache_DISCRETIZE_TEW___['qn_BF_dy_x1']
            qn_BF_dy_y1 = self.___cache_DISCRETIZE_TEW___['qn_BF_dy_y1']
            quad_weights = self.___cache_DISCRETIZE_TEW___['quad_weights']

        #------- check func and get func --------------------------------------------------------------------
        if target == 'func':
            assert self.func.body is not None, f"No func.body!"
            TEW_func = self.func.body
        else:
            raise NotImplementedError(f"1Trace = discretize_VectorField_standard of ftype trace-element-wise "
                                      f"not applicable for target={target}.")

        # --- Now, we do a check whether trace-elements in TEW_func are local -------------------------------
        for T in TEW_func:
            assert T in self.mesh.trace.elements, f"trace-element #{T} is not a local trace-element."

        # dispatch the lens for integration ------------------------------------------------
        lens_NS_dy, lens_NS_dz, lens_WE_dx, lens_WE_dz, lens_BF_dx, lens_BF_dy = LENS
        # initialize the Trace-Element-Wise-local-cochain dict ------------------------------------
        local_TEW = dict()
        Zo = [0,]
        # Go through all valid local trace-elements in the function ---------------------------
        for T in TEW_func:
            te = self.mesh.trace.elements[T]
            # ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side

            #--------------------------------- NS sides ---------------------------------------
            if ele_side in 'NS':
                J = te.coordinate_transformation.Jacobian_matrix(qn_NS_dy_y, qn_NS_dy_z)
                J = (J[0][0], J[1][0], J[2][0]) # dy of (dy, dz)
                u, v, w = list(), list(), list()
                for dy, dz in zip(qn_NS_dy_y1, qn_NS_dy_z1):
                    ___, uvw = TEW_func[T](Zo, dy, dz)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,:,0]
                v = np.array(v)[:,:,0]
                w = np.array(w)[:,:,0]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[1]
                C = lens_NS_dy
                cochain_dy = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                J = te.coordinate_transformation.Jacobian_matrix(qn_NS_dz_y, qn_NS_dz_z)
                J = (J[0][1], J[1][1], J[2][1]) # dz of (dy, dz)
                u, v, w = list(), list(), list()
                for dy, dz in zip(qn_NS_dz_y1, qn_NS_dz_z1):
                    ___, uvw = TEW_func[T](Zo, dy, dz)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,0,:]
                v = np.array(v)[:,0,:]
                w = np.array(w)[:,0,:]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[2]
                C = lens_NS_dz
                cochain_dz = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dy, cochain_dz])

            #--------------------------------- WE sides ---------------------------------------
            elif ele_side in 'WE':
                J = te.coordinate_transformation.Jacobian_matrix(qn_WE_dx_x, qn_WE_dx_z)
                J = (J[0][1], J[1][1], J[2][1]) # dx of (dz, dx)
                u, v, w = list(), list(), list()
                for dx, dz in zip(qn_WE_dx_x1, qn_WE_dx_z1):
                    ___, uvw = TEW_func[T](dx, Zo, dz)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,:,0]
                v = np.array(v)[:,:,0]
                w = np.array(w)[:,:,0]

                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[0]
                C = lens_WE_dx
                cochain_dx = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                J = te.coordinate_transformation.Jacobian_matrix(qn_WE_dz_x, qn_WE_dz_z)
                J = (J[0][0], J[1][0], J[2][0]) # dz of (dz, dx)
                u, v, w = list(), list(), list()
                for dx, dz in zip(qn_WE_dz_x1, qn_WE_dz_z1):
                    ___, uvw = TEW_func[T](dx, Zo, dz)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,0,:]
                v = np.array(v)[:,0,:]
                w = np.array(w)[:,0,:]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[2]
                C = lens_WE_dz
                cochain_dz = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dx, cochain_dz])

            #--------------------------------- BF sides ---------------------------------------
            elif ele_side in 'BF':

                J = te.coordinate_transformation.Jacobian_matrix(qn_BF_dx_x, qn_BF_dx_y)
                J = (J[0][0], J[1][0], J[2][0]) # dx of (dx, dy)
                u, v, w = list(), list(), list()
                for dx, dy in zip(qn_BF_dx_x1, qn_BF_dx_y1):
                    ___, uvw = TEW_func[T](dx, dy, Zo)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,:,0]
                v = np.array(v)[:,:,0]
                w = np.array(w)[:,:,0]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[0]
                C = lens_BF_dx
                cochain_dx = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                J = te.coordinate_transformation.Jacobian_matrix(qn_BF_dy_x, qn_BF_dy_y)
                J = (J[0][1], J[1][1], J[2][1]) # dy of (dx, dy)
                u, v, w = list(), list(), list()
                for dx, dy in zip(qn_BF_dy_x1, qn_BF_dy_y1):
                    ___, uvw = TEW_func[T](dx, dy, Zo)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,0,:]
                v = np.array(v)[:,0,:]
                w = np.array(w)[:,0,:]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[1]
                C = lens_BF_dy
                cochain_dy = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dx, cochain_dy])

            else:
                raise Exception()

            if not self.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[T] = te_primal_local

        if update_cochain: self.cochain.local_TEW = local_TEW
        # 'locally full local TEW cochain': provide cochain.local_TEW and for all dofs on the trace element.
        return 'locally full local TEW cochain', local_TEW



    def reconstruct(self, xi, eta, sigma, ravel=False, i=None):
        """
        Do the reconstruction.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param sigma: A 1d iterable object of floats between -1 and 1.
        :param bool ravel: (`default`:``False``) If we return 1d data?
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :type sigma: list, tuple, numpy.ndarray
        :param i: (`default`:``None``) Do the reconstruction for these
            trace elements. if it is ``None``, then do it for all trace
            elements.
        :type i: int, None, list, tuple
        """
        if i is None:
            indices = self.mesh.trace.elements._elements_.keys()
        else:
            if not isinstance(i, (list, tuple)):
                indices = [i,]
            else:
                indices = i

        px, py, pz = self.p
        # N, S trace element
        num_NS_dy = py*(pz+1)
        num_NS_dz = pz*(py+1)
        # W, E trace element
        num_WE_dx = px*(pz+1)
        num_WE_dz = pz*(px+1)
        # B, F trace element
        num_BF_dx = px*(py+1)
        num_BF_dy = py*(px+1)

        xietasigma, pb = self.DO.evaluate_basis_at_meshgrid(xi, eta, sigma)
        ii, jj, kk = np.size(xi), np.size(eta), np.size(sigma)
        xyz = dict()
        v = dict()

        for key in indices:
            if key in self.mesh.trace.elements:
                te = self.mesh.trace.elements[key]
                side = te.CHARACTERISTIC_side
                ele = te.CHARACTERISTIC_element
                xyz_i = te.coordinate_transformation.mapping(*xietasigma[side], from_element=ele, side=side)
                iJ = te.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma[side])
                prime_cochain = self.cochain.local_TEW[key]

                b0, b1 = pb[side]
                if side in 'NS':
                    assert num_NS_dy + num_NS_dz == len(prime_cochain), "A trivial check!"
                    c0 = prime_cochain[:num_NS_dy]
                    c1 = prime_cochain[-num_NS_dz:]
                    iJ00, iJ01 = iJ[0][0], iJ[1][0]
                    iJ10, iJ11 = iJ[0][1], iJ[1][1]
                    iJ20, iJ21 = iJ[0][2], iJ[1][2]
                elif side in 'WE':
                    assert num_WE_dx + num_WE_dz == len(prime_cochain), "A trivial check!"
                    c0 = prime_cochain[:num_WE_dx]
                    c1 = prime_cochain[-num_WE_dz:]
                    iJ01, iJ00 = iJ[0][0], iJ[1][0]
                    iJ11, iJ10 = iJ[0][1], iJ[1][1]
                    iJ21, iJ20 = iJ[0][2], iJ[1][2]
                elif side in 'BF':
                    assert num_BF_dx + num_BF_dy == len(prime_cochain), "A trivial check!"
                    c0 = prime_cochain[:num_BF_dx]
                    c1 = prime_cochain[-num_BF_dy:]
                    iJ00, iJ01 = iJ[0][0], iJ[1][0]
                    iJ10, iJ11 = iJ[0][1], iJ[1][1]
                    iJ20, iJ21 = iJ[0][2], iJ[1][2]
                else:
                    raise Exception()

                v0 = np.einsum('ij, i -> j', b0, c0, optimize='greedy')
                v1 = np.einsum('ij, i -> j', b1, c1, optimize='greedy')

                v_x = v0*iJ00 + v1*iJ01
                v_y = v0*iJ10 + v1*iJ11
                v_z = v0*iJ20 + v1*iJ21

                if ravel:
                    xyz[key] = xyz_i
                    v[key] = [v_x, v_y, v_z]
                else:
                    if side in 'NS':
                        xyz[key] = [xyz_i[m].reshape(jj, kk, order='F') for m in range(3)]
                        v[key] = [v_x.reshape((jj, kk), order='F'), v_y.reshape((jj, kk), order='F'), v_z.reshape((jj, kk), order='F')]
                    elif side in 'WE':
                        xyz[key] = [xyz_i[m].reshape(ii, kk, order='F') for m in range(3)]
                        v[key] = [v_x.reshape((ii, kk), order='F'), v_y.reshape((ii, kk), order='F'), v_z.reshape((ii, kk), order='F')]
                    elif side in 'BF':
                        xyz[key] = [xyz_i[m].reshape(ii, jj, order='F') for m in range(3)]
                        v[key] = [v_x.reshape((ii, jj), order='F'), v_y.reshape((ii, jj), order='F'), v_z.reshape((ii, jj), order='F')]
                    else:
                        raise Exception

        return xyz, v



    def ___PRIVATE_generate_TEW_mass_matrices___(self):
        """Generate the trace-element-wise mass matrices stored in a dict whose keys are trace-element numbers
        and values are the mass matrices in the corresponding trace-elements.
        """
        p = [self.dqp[i]+2 for i in range(self.ndim)] # +2 for safety, the mass matrices of standard forms use dqp
        quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad

        qw = dict()
        qw['NS'] = np.kron(quad_weights[2], quad_weights[1])
        qw['WE'] = np.kron(quad_weights[2], quad_weights[0])
        qw['BF'] = np.kron(quad_weights[1], quad_weights[0])

        xietasigma, pb = self.DO.evaluate_basis_at_meshgrid(*quad_nodes)

        local_cache = dict()

        MD = dict()
        for i in self.mesh.trace.elements:
            te = self.mesh.trace.elements[i]
            side = te.CHARACTERISTIC_side
            mark = te.type_wrt_metric.mark

            if isinstance(mark, str) and mark in local_cache: # not an id (chaotic) mark, could be cached.
                MD[i] = local_cache[mark]
            else:
                g = te.coordinate_transformation.metric(*xietasigma[side])
                iG = te.coordinate_transformation.inverse_metric_matrix(*xietasigma[side])
                b0, b1 = pb[side]
                if side in 'NS':
                    QW = qw['NS']
                elif side in 'WE':
                    QW = qw['WE']
                elif side in 'BF':
                    QW = qw['BF']
                else:
                    raise Exception()

                if side in 'NSBF':

                    M00 = np.einsum('im, jm, m -> ij', b0, b0, np.sqrt(g) * iG[0][0] * QW, optimize='greedy')
                    M11 = np.einsum('im, jm, m -> ij', b1, b1, np.sqrt(g) * iG[1][1] * QW, optimize='greedy')
                    if isinstance(mark, str) and mark[:5] == 'Orth.':
                        M01 = np.zeros((M00.shape[0], M11.shape[1]))
                    else:
                        M01 = np.einsum('im, jm, m -> ij', b0, b1, np.sqrt(g) * iG[0][1] * QW, optimize='greedy')

                else: # WE sides
                    M00 = np.einsum('im, jm, m -> ij', b0, b0, np.sqrt(g) * iG[1][1] * QW, optimize='greedy')
                    M11 = np.einsum('im, jm, m -> ij', b1, b1, np.sqrt(g) * iG[0][0] * QW, optimize='greedy')
                    if isinstance(mark, str) and mark[:5] == 'Orth.':
                        M01 = np.zeros((M00.shape[0], M11.shape[1]))
                    else:
                        M01 = np.einsum('im, jm, m -> ij', b0, b1, np.sqrt(g) * iG[1][0] * QW, optimize='greedy')

                M10 = M01.T
                M = spspa.csc_matrix(np.bmat([(M00, M01), (M10, M11)]))

                # print(iG[1][1])

                if isinstance(mark, str): local_cache[mark] = M

                MD[i] = M

        return MD












if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\form\trace\_1_trace.py

    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.1)([2,3,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',4), ('Lobatto',6)])
    FC = FormCaller(mesh, space)

    def u(t, x, y, z): return t + np.sin(2*np.pi*x) * np.cos(np.pi*y) * np.cos(2*np.pi*z)
    def v(t, x, y, z): return t + np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(np.pi*z)
    def w(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z)
    # def u(t, x, y, z): return 0.3 + 0*x
    # def v(t, x, y, z): return 0.75 + 0*x
    # def w(t, x, y, z): return 0.25 + 0*x
    V = FC('vector', (u,v,w))

    t1 = FC('1-t')

    t1.TW.func.DO.set_func_body_as(V)
    t1.TW.current_time = 0
    t1.TW.DO.push_all_to_instant()

    t1.discretize()

    #
    # xi, et, sg = np.linspace(-1,1,9), np.linspace(-1,1,10), np.linspace(-1,1,11)
    #
    # xyz, v = t1.reconstruct(xi, et, sg)
    # print(v.keys())

    # t1.visualize()

    tM1 = t1.matrices.mass





