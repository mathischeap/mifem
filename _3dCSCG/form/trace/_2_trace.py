# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config import *
from root.mifem import read
from scipy.interpolate import NearestNDInterpolator
from SCREWS.quadrature import Quadrature
from _3dCSCG.form.trace.base.main import _3dCSCG_Standard_Trace


class _2Trace(_3dCSCG_Standard_Trace):
    """
    Trace 2-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-2-trace-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_trace_2form')
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self.___cache_DISCRETIZE_STANDARD___ = None
        super().RESET_cache()




    def ___TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3


        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 2-trace FUNC cannot accommodate _3dCSCG_ScalarField of ftype {func_body.ftype}."
        elif func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 2-trace FUNC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise NotImplementedError(
                f"3d CSCG 2-trace form FUNC cannot accommodate {func_body}.")

    def ___TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3


        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard', 'boundary-wise'), \
                f"3dCSCG 2-trace BC cannot accommodate _3dCSCG_ScalarField of ftype {func_body.ftype}."
        elif func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard', 'boundary-wise'), \
                f"3dCSCG 2-trace BC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise NotImplementedError(
                f"3d CSCG 2-trace form BC cannot accommodate {func_body}.")




    def discretize(self, update_cochain=True, target='func', **kwargs):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param target:
        :param kwargs: Keywords arguments to be passed to particular discretization schemes.
        :return: The cochain corresponding to the particular discretization scheme.
        """
        if target == 'func':
            if self.TW.func.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_ScalarField_of_ftype_standard___(
                        update_cochain=update_cochain, **kwargs)
                else:
                    raise Exception(f'3dCSCG 2-trace can not (target func) discretize _3dCSCG_ScalarField of ftype {self.func.ftype}.')

            elif self.TW.func.body.__class__.__name__ == '_3dCSCG_VectorField':
                if self.func.ftype == 'standard': # we will discretize the norm component of the vector.
                    return self.___PRIVATE_discretize_the_flux_of_a_VectorField_of_ftype_standard___(
                        update_cochain=update_cochain, **kwargs)
                else:
                    raise Exception(f'3dCSCG 2-trace can not (target func) discretize _3dCSCG_VectorField of ftype {self.func.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 2-trace can not (target func) discretize {self.TW.func.body.__class__}.')

        elif target == 'BC': # We target at the BC, so we do not update the cochain!

            if self.TW.BC.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if self.BC.ftype == 'standard':
                    return self.___PRIVATE_discretize_ScalarField_of_ftype_standard___(
                        update_cochain=False, target='BC', **kwargs)
                elif self.BC.ftype == 'boundary-wise':
                    return self.___PRIVATE_discretize_ScalarField_of_ftype_boundary_wise___(
                        **kwargs) # must be False update_cochain and 'BC' target.
                else:
                    raise Exception(f'3dCSCG 2-trace can not (target BC) discretize _3dCSCG_ScalarField of ftype {self.BC.ftype}.')

            elif self.TW.BC.body.__class__.__name__ == '_3dCSCG_VectorField':
                if self.BC.ftype == 'standard': # we will discretize the norm flux of the vector.
                    return self.___PRIVATE_discretize_the_flux_of_a_VectorField_of_ftype_standard___(
                        update_cochain=False, target='BC', **kwargs)
                elif self.BC.ftype == 'boundary-wise': # we will discretize the norm flux of the vector.
                    return self.___PRIVATE_discretize_the_flux_of_a_VectorField_of_ftype_boundary_wise___(
                        **kwargs) # must be False update_cochain and 'BC' target.
                else:
                    raise Exception(f'3dCSCG 2-trace can not (target BC) discretize _3dCSCG_VectorField of ftype {self.BC.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 2-trace can not (target BC) discretize {self.TW.BC.body.__class__}.')

        else:
            raise NotImplementedError(f"target={target} not implemented "
                                      f"for 3d CSCG 2-trace form discretization.")

    def ___PRIVATE_discretize_ScalarField_of_ftype_standard___(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the standard scalar field to all trace elements.

        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.
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
            # NS sides
            qn_NS_y = []
            qn_NS_z = []
            for k in range(self.p[2]):
                for j in range(self.p[1]):
                    qn_NS_y.append(qnodes[1][j][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_NS_z.append(qnodes[2][k][np.newaxis,:].repeat(p[1]+1, axis=0))
            qn_NS_y, qn_NS_z = np.array(qn_NS_y), np.array(qn_NS_z)
            area_NS = np.kron(lens[2], lens[1]) * 0.25
            # WE sides
            qn_WE_x = []
            qn_WE_z = []
            for k in range(self.p[2]):
                for i in range(self.p[0]):
                    qn_WE_x.append(qnodes[0][i][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_WE_z.append(qnodes[2][k][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_WE_x, qn_WE_z = np.array(qn_WE_x), np.array(qn_WE_z)
            area_WE = np.kron(lens[2], lens[0]) * 0.25
            # BF sides
            qn_BF_x = []
            qn_BF_y = []
            for j in range(self.p[1]):
                for i in range(self.p[0]):
                    qn_BF_x.append(qnodes[0][i][:,np.newaxis].repeat(p[1]+1, axis=1))
                    qn_BF_y.append(qnodes[1][j][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_BF_x, qn_BF_y = np.array(qn_BF_x), np.array(qn_BF_y)
            area_BF = np.kron(lens[1], lens[0]) * 0.25

            cd = dict()
            cd['quadDegree'] = quad_degree
            cd['qn_NS_y'] = qn_NS_y
            cd['qn_NS_z'] = qn_NS_z
            cd['area_NS'] = area_NS
            cd['qn_WE_x'] = qn_WE_x
            cd['qn_WE_z'] = qn_WE_z
            cd['area_WE'] = area_WE
            cd['qn_BF_x'] = qn_BF_x
            cd['qn_BF_y'] = qn_BF_y
            cd['area_BF'] = area_BF
            cd['quad_weights'] = quad_weights
            self.___cache_DISCRETIZE_STANDARD___ = cd
        else:
            qn_NS_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_y']
            qn_NS_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_z']
            area_NS = self.___cache_DISCRETIZE_STANDARD___['area_NS']
            qn_WE_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_x']
            qn_WE_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_z']
            area_WE = self.___cache_DISCRETIZE_STANDARD___['area_WE']
            qn_BF_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_x']
            qn_BF_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_y']
            area_BF = self.___cache_DISCRETIZE_STANDARD___['area_BF']
            quad_weights = self.___cache_DISCRETIZE_STANDARD___['quad_weights']

        if target == 'func':
            assert self.func.body is not None, f"No func.body!"
            _lf_ = self.func.body[0]
        elif target == 'BC':
            assert self.BC.body is not None, f"No BC.body!"
            _lf_ = self.BC.body[0]
        else:
            raise NotImplementedError(f"Not applicable for target={target}.")

        local_TEW = dict()
        for key in self.mesh.trace.elements:
            te = self.mesh.trace.elements[key]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side
            if ele_side in 'NS':
                qn0, qn1 = qn_NS_y, qn_NS_z
                qw0, qw1 = quad_weights[1], quad_weights[2]
                area = area_NS
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
            elif ele_side in 'WE':
                qn0, qn1 = qn_WE_x, qn_WE_z
                qw0, qw1 = quad_weights[0], quad_weights[2]
                area = area_WE
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
            elif ele_side in 'BF':
                qn0, qn1 = qn_BF_x, qn_BF_y
                qw0, qw1 = quad_weights[0], quad_weights[1]
                area = area_BF
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
            else:
                raise Exception()

            f = _lf_(x, y, z)
            sqrt_g = np.sqrt(g)
            te_primal_local = np.einsum('mij, i, j, m -> m', f * sqrt_g,
                                        qw0, qw1, area,
                                        optimize='greedy')

            if not self.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[key] = te_primal_local

        if update_cochain: self.cochain.local_TEW = local_TEW
        # 'locally full local TEW cochain': provide cochain.local_TEW and for all dofs on the trace element.
        return 'locally full local TEW cochain', local_TEW

    def ___PRIVATE_discretize_ScalarField_of_ftype_boundary_wise___(self, quad_degree=None):
        """"""

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
            # NS sides
            qn_NS_y = []
            qn_NS_z = []
            for k in range(self.p[2]):
                for j in range(self.p[1]):
                    qn_NS_y.append(qnodes[1][j][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_NS_z.append(qnodes[2][k][np.newaxis,:].repeat(p[1]+1, axis=0))
            qn_NS_y, qn_NS_z = np.array(qn_NS_y), np.array(qn_NS_z)
            area_NS = np.kron(lens[2], lens[1]) * 0.25
            # WE sides
            qn_WE_x = []
            qn_WE_z = []
            for k in range(self.p[2]):
                for i in range(self.p[0]):
                    qn_WE_x.append(qnodes[0][i][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_WE_z.append(qnodes[2][k][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_WE_x, qn_WE_z = np.array(qn_WE_x), np.array(qn_WE_z)
            area_WE = np.kron(lens[2], lens[0]) * 0.25
            # BF sides
            qn_BF_x = []
            qn_BF_y = []
            for j in range(self.p[1]):
                for i in range(self.p[0]):
                    qn_BF_x.append(qnodes[0][i][:,np.newaxis].repeat(p[1]+1, axis=1))
                    qn_BF_y.append(qnodes[1][j][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_BF_x, qn_BF_y = np.array(qn_BF_x), np.array(qn_BF_y)
            area_BF = np.kron(lens[1], lens[0]) * 0.25

            cd = dict()
            cd['quadDegree'] = quad_degree
            cd['qn_NS_y'] = qn_NS_y
            cd['qn_NS_z'] = qn_NS_z
            cd['area_NS'] = area_NS
            cd['qn_WE_x'] = qn_WE_x
            cd['qn_WE_z'] = qn_WE_z
            cd['area_WE'] = area_WE
            cd['qn_BF_x'] = qn_BF_x
            cd['qn_BF_y'] = qn_BF_y
            cd['area_BF'] = area_BF
            cd['quad_weights'] = quad_weights
            self.___cache_DISCRETIZE_STANDARD___ = cd
        else:
            qn_NS_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_y']
            qn_NS_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_z']
            area_NS = self.___cache_DISCRETIZE_STANDARD___['area_NS']
            qn_WE_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_x']
            qn_WE_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_z']
            area_WE = self.___cache_DISCRETIZE_STANDARD___['area_WE']
            qn_BF_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_x']
            qn_BF_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_y']
            area_BF = self.___cache_DISCRETIZE_STANDARD___['area_BF']
            quad_weights = self.___cache_DISCRETIZE_STANDARD___['quad_weights']

        assert self.BC.body is not None, f"No BC.body!"
        FUNC = self.BC.body
        RANGE_trace_elements = self.mesh.boundaries.RANGE_trace_elements
        local_TEW = dict()
        for bn in FUNC:
            func_bn = FUNC[bn]
            trace_elements = RANGE_trace_elements[bn]
            _lf_ = func_bn[0]
            for i in trace_elements:
                te = self.mesh.trace.elements[i]
                ele = te.CHARACTERISTIC_element
                ele_side = te.CHARACTERISTIC_side
                if ele_side in 'NS':
                    qn0, qn1 = qn_NS_y, qn_NS_z
                    qw0, qw1 = quad_weights[1], quad_weights[2]
                    area = area_NS
                    x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                    g = te.coordinate_transformation.metric(qn0, qn1)
                elif ele_side in 'WE':
                    qn0, qn1 = qn_WE_x, qn_WE_z
                    qw0, qw1 = quad_weights[0], quad_weights[2]
                    area = area_WE
                    x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                    g = te.coordinate_transformation.metric(qn0, qn1)
                elif ele_side in 'BF':
                    qn0, qn1 = qn_BF_x, qn_BF_y
                    qw0, qw1 = quad_weights[0], quad_weights[1]
                    area = area_BF
                    x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                    g = te.coordinate_transformation.metric(qn0, qn1)
                else:
                    raise Exception()

                f = _lf_(x, y, z)
                sqrt_g = np.sqrt(g)
                te_primal_local = np.einsum('mij, i, j, m -> m', f * sqrt_g,
                                            qw0, qw1, area,
                                            optimize='greedy')
                if not self.space.IS_Kronecker: raise NotImplementedError()
                assert i not in local_TEW, f"Trace element #{i} can only appear once (be on one mesh boundary)."
                local_TEW[i] = te_primal_local

        # 'locally full local TEW cochain': provide cochain.local_TEW and for all dofs on the trace element.
        return 'locally full local TEW cochain', local_TEW


    def ___PRIVATE_discretize_the_flux_of_a_VectorField_of_ftype_standard___(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the norm component of the vector to the trace 2-form.

        As a trace 2-form can only accommodate scalar, when given a
        vector, for example denoted by vec(u), we actually mean we get
        vec(u) · vec(n) where vec(n) represent the positive normal
        direction of the surface (trace element). The positive normal
        direction means x+, y+ or z+ direction, or more accurately,
        xi+, eta+ or sigma+ direction (in the mesh element setting)
        before mapping.


        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.
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
            # NS sides
            qn_NS_y = []
            qn_NS_z = []
            for k in range(self.p[2]):
                for j in range(self.p[1]):
                    qn_NS_y.append(qnodes[1][j][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_NS_z.append(qnodes[2][k][np.newaxis,:].repeat(p[1]+1, axis=0))
            qn_NS_y, qn_NS_z = np.array(qn_NS_y), np.array(qn_NS_z)
            area_NS = np.kron(lens[2], lens[1]) * 0.25
            # WE sides
            qn_WE_x = []
            qn_WE_z = []
            for k in range(self.p[2]):
                for i in range(self.p[0]):
                    qn_WE_x.append(qnodes[0][i][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_WE_z.append(qnodes[2][k][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_WE_x, qn_WE_z = np.array(qn_WE_x), np.array(qn_WE_z)
            area_WE = np.kron(lens[2], lens[0]) * 0.25
            # BF sides
            qn_BF_x = []
            qn_BF_y = []
            for j in range(self.p[1]):
                for i in range(self.p[0]):
                    qn_BF_x.append(qnodes[0][i][:,np.newaxis].repeat(p[1]+1, axis=1))
                    qn_BF_y.append(qnodes[1][j][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_BF_x, qn_BF_y = np.array(qn_BF_x), np.array(qn_BF_y)
            area_BF = np.kron(lens[1], lens[0]) * 0.25

            cd = dict()
            cd['quadDegree'] = quad_degree
            cd['qn_NS_y'] = qn_NS_y
            cd['qn_NS_z'] = qn_NS_z
            cd['area_NS'] = area_NS
            cd['qn_WE_x'] = qn_WE_x
            cd['qn_WE_z'] = qn_WE_z
            cd['area_WE'] = area_WE
            cd['qn_BF_x'] = qn_BF_x
            cd['qn_BF_y'] = qn_BF_y
            cd['area_BF'] = area_BF
            cd['quad_weights'] = quad_weights
            self.___cache_DISCRETIZE_STANDARD___ = cd
        else:
            qn_NS_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_y']
            qn_NS_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_z']
            area_NS = self.___cache_DISCRETIZE_STANDARD___['area_NS']
            qn_WE_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_x']
            qn_WE_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_z']
            area_WE = self.___cache_DISCRETIZE_STANDARD___['area_WE']
            qn_BF_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_x']
            qn_BF_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_y']
            area_BF = self.___cache_DISCRETIZE_STANDARD___['area_BF']
            quad_weights = self.___cache_DISCRETIZE_STANDARD___['quad_weights']

        if target == 'func':
            assert self.func.body is not None, f"No func.body!"
            _lf0_, _lf1_, _lf2_ = self.func.body
        elif target == 'BC':
            assert self.BC.body is not None, f"No BC.body!"
            _lf0_, _lf1_, _lf2_ = self.BC.body
        else:
            raise NotImplementedError(f"Not applicable for target={target}.")

        local_TEW = dict()
        for key in self.mesh.trace.elements:
            te = self.mesh.trace.elements[key]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side
            if ele_side in 'NS':
                qn0, qn1 = qn_NS_y, qn_NS_z
                qw0, qw1 = quad_weights[1], quad_weights[2]
                area = area_NS
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
                ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
            elif ele_side in 'WE':
                qn0, qn1 = qn_WE_x, qn_WE_z
                qw0, qw1 = quad_weights[0], quad_weights[2]
                area = area_WE
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
                ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
            elif ele_side in 'BF':
                qn0, qn1 = qn_BF_x, qn_BF_y
                qw0, qw1 = quad_weights[0], quad_weights[1]
                area = area_BF
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
                ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
            else:
                raise Exception()

            f0, f1, f2 = _lf0_(x, y, z), _lf1_(x, y, z), _lf2_(x, y, z)
            f = f0 * ouv[0] + f1 * ouv[1] + f2 * ouv[2]
            sqrt_g = np.sqrt(g)
            te_primal_local = np.einsum('mij, i, j, m -> m', f * sqrt_g,
                                        qw0, qw1, area,
                                        optimize='greedy')

            if not self.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[key] = te_primal_local


        if update_cochain: self.cochain.local_TEW = local_TEW
        # 'locally full local TEW cochain': provide cochain.local_TEW and for all dofs on the trace element.
        return 'locally full local TEW cochain', local_TEW

    def ___PRIVATE_discretize_the_flux_of_a_VectorField_of_ftype_boundary_wise___(self, quad_degree=None):
        """We will discretize the norm component of the vector to the trace 2-form.

        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.
        """

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
            # NS sides
            qn_NS_y = []
            qn_NS_z = []
            for k in range(self.p[2]):
                for j in range(self.p[1]):
                    qn_NS_y.append(qnodes[1][j][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_NS_z.append(qnodes[2][k][np.newaxis,:].repeat(p[1]+1, axis=0))
            qn_NS_y, qn_NS_z = np.array(qn_NS_y), np.array(qn_NS_z)
            area_NS = np.kron(lens[2], lens[1]) * 0.25
            # WE sides
            qn_WE_x = []
            qn_WE_z = []
            for k in range(self.p[2]):
                for i in range(self.p[0]):
                    qn_WE_x.append(qnodes[0][i][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_WE_z.append(qnodes[2][k][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_WE_x, qn_WE_z = np.array(qn_WE_x), np.array(qn_WE_z)
            area_WE = np.kron(lens[2], lens[0]) * 0.25
            # BF sides
            qn_BF_x = []
            qn_BF_y = []
            for j in range(self.p[1]):
                for i in range(self.p[0]):
                    qn_BF_x.append(qnodes[0][i][:,np.newaxis].repeat(p[1]+1, axis=1))
                    qn_BF_y.append(qnodes[1][j][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_BF_x, qn_BF_y = np.array(qn_BF_x), np.array(qn_BF_y)
            area_BF = np.kron(lens[1], lens[0]) * 0.25

            cd = dict()
            cd['quadDegree'] = quad_degree
            cd['qn_NS_y'] = qn_NS_y
            cd['qn_NS_z'] = qn_NS_z
            cd['area_NS'] = area_NS
            cd['qn_WE_x'] = qn_WE_x
            cd['qn_WE_z'] = qn_WE_z
            cd['area_WE'] = area_WE
            cd['qn_BF_x'] = qn_BF_x
            cd['qn_BF_y'] = qn_BF_y
            cd['area_BF'] = area_BF
            cd['quad_weights'] = quad_weights
            self.___cache_DISCRETIZE_STANDARD___ = cd
        else:
            qn_NS_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_y']
            qn_NS_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_z']
            area_NS = self.___cache_DISCRETIZE_STANDARD___['area_NS']
            qn_WE_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_x']
            qn_WE_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_z']
            area_WE = self.___cache_DISCRETIZE_STANDARD___['area_WE']
            qn_BF_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_x']
            qn_BF_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_y']
            area_BF = self.___cache_DISCRETIZE_STANDARD___['area_BF']
            quad_weights = self.___cache_DISCRETIZE_STANDARD___['quad_weights']


        assert self.BC.body is not None, f"No BC.body!"
        FUNC = self.BC.body

        RANGE_trace_elements = self.mesh.boundaries.RANGE_trace_elements

        local_TEW = dict()
        for bn in FUNC:
            func_bn = FUNC[bn]
            trace_elements = RANGE_trace_elements[bn]
            _lf0_, _lf1_, _lf2_ = func_bn
            for i in trace_elements:
                te = self.mesh.trace.elements[i]
                ele = te.CHARACTERISTIC_element
                ele_side = te.CHARACTERISTIC_side
                if ele_side in 'NS':
                    qn0, qn1 = qn_NS_y, qn_NS_z
                    qw0, qw1 = quad_weights[1], quad_weights[2]
                    area = area_NS
                    x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                    g = te.coordinate_transformation.metric(qn0, qn1)
                    ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
                elif ele_side in 'WE':
                    qn0, qn1 = qn_WE_x, qn_WE_z
                    qw0, qw1 = quad_weights[0], quad_weights[2]
                    area = area_WE
                    x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                    g = te.coordinate_transformation.metric(qn0, qn1)
                    ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
                elif ele_side in 'BF':
                    qn0, qn1 = qn_BF_x, qn_BF_y
                    qw0, qw1 = quad_weights[0], quad_weights[1]
                    area = area_BF
                    x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                    g = te.coordinate_transformation.metric(qn0, qn1)
                    ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
                else:
                    raise Exception()

                f0, f1, f2 = _lf0_(x, y, z), _lf1_(x, y, z), _lf2_(x, y, z)
                f = f0 * ouv[0] + f1 * ouv[1] + f2 * ouv[2]
                sqrt_g = np.sqrt(g)
                te_primal_local = np.einsum('mij, i, j, m -> m',
                                            f * sqrt_g,
                                            qw0, qw1, area,
                                            optimize='greedy')
                if not self.space.IS_Kronecker: raise NotImplementedError()
                assert i not in local_TEW, f"Trace element #{i} can only appear once (be on one mesh boundary)."
                local_TEW[i] = te_primal_local

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
                g = te.coordinate_transformation.metric(*xietasigma[side])
                prime_cochain = self.cochain.local_TEW[key]
                if side in 'NS':
                    vi = np.einsum('i, j, ij -> j', prime_cochain, 1 / np.sqrt(g),
                                   pb[side][0], optimize='greedy')
                elif side in 'WE':
                    vi = np.einsum('i, j, ij -> j', prime_cochain, 1 / np.sqrt(g),
                                   pb[side][0], optimize='greedy')
                elif side in 'BF':
                    vi = np.einsum('i, j, ij -> j', prime_cochain, 1 / np.sqrt(g),
                                   pb[side][0], optimize='greedy')
                else:
                    raise Exception()
                if ravel:
                    xyz[key] = xyz_i
                    v[key] = [vi,]
                else:
                    if side in 'NS':
                        xyz[key] = [xyz_i[m].reshape(jj, kk, order='F') for m in range(3)]
                        v[key] = [vi.reshape((jj, kk), order='F'),]
                    elif side in 'WE':
                        xyz[key] = [xyz_i[m].reshape(ii, kk, order='F') for m in range(3)]
                        v[key] = [vi.reshape((ii, kk), order='F'),]
                    elif side in 'BF':
                        xyz[key] = [xyz_i[m].reshape(ii, jj, order='F') for m in range(3)]
                        v[key] = [vi.reshape((ii, jj), order='F'),]
                    else:
                        raise Exception
        return xyz, v


    def ___PRIVATE_generate_TEW_mass_matrices___(self):
        """Generate the trace-element-wise mass matrices."""
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
                b = pb[side][0]

                if side in 'NS':
                    M = np.einsum('im, jm, m -> ij', b, b, np.reciprocal(np.sqrt(g)) * qw['NS'], optimize='greedy')
                elif side in 'WE':
                    M = np.einsum('im, jm, m -> ij', b, b, np.reciprocal(np.sqrt(g)) * qw['WE'], optimize='greedy')
                elif side in 'BF':
                    M = np.einsum('im, jm, m -> ij', b, b, np.reciprocal(np.sqrt(g)) * qw['BF'], optimize='greedy')
                else:
                    raise Exception()

                if isinstance(mark, str): local_cache[mark] = M

                MD[i] = M

        return MD


    def ___DO_resemble___(self, obj_or_filename):
        """

        :param obj_or_filename:
        :return:
        """
        if isinstance(obj_or_filename, str):
            ot = read(obj_or_filename)
        else:
            ot = obj_or_filename

        assert ot.mesh.domain == self.mesh.domain, "domain must be same."
        assert self.__class__.__name__ == ot.__class__.__name__
        assert self.mesh.__class__.__name__ == ot.mesh.__class__.__name__

        bp = int(np.ceil((20000 / self.mesh.elements.GLOBAL_num) ** (1/3)))
        p = [bp + self.p[i] for i in range(3)]
        gap = [1 / (p[i]+1) for i in range(3)]
        r = np.linspace(-1 + gap[0], 1 - gap[0], p[0])
        s = np.linspace(-1 + gap[1], 1 - gap[1], p[1])
        t = np.linspace(-1 + gap[2], 1 - gap[2], p[2])
        xyz, V = ot.reconstruct(r, s, t, ravel=True)

        xyz = cOmm.gather(xyz, root=mAster_rank)
        V = cOmm.gather(V, root=mAster_rank)

        tep = dict()
        for i in ot.mesh.trace.elements:
            TEi = ot.mesh.trace.elements[i]
            tep[i] = TEi.CHARACTERISTIC_side
        tep = cOmm.gather(tep, root=mAster_rank)

        if rAnk == mAster_rank:
            XYZ, VVV, TEP = dict(), dict(), dict()
            for i in range(len(xyz)):
                XYZ.update(xyz[i])
                VVV.update(V[i])
                TEP.update(tep[i])
            del xyz, V, tep
            V_x, x_x, x_y, x_z, V_y, y_x, y_y, y_z, V_z, z_x, z_y, z_z = \
                [np.array([]) for _ in range(12)]
            I_func = dict()
            for i in range(ot.mesh.trace.elements.GLOBAL_num):
                ele_side = TEP[i]
                if ele_side in 'NS':
                    x_x = np.append(x_x, XYZ[i][0])
                    x_y = np.append(x_y, XYZ[i][1])
                    x_z = np.append(x_z, XYZ[i][2])
                    V_x = np.append(V_x, VVV[i][0])
                elif ele_side in 'WE':
                    y_x = np.append(y_x, XYZ[i][0])
                    y_y = np.append(y_y, XYZ[i][1])
                    y_z = np.append(y_z, XYZ[i][2])
                    V_y = np.append(V_y, VVV[i][0])
                elif ele_side in 'BF':
                    z_x = np.append(z_x, XYZ[i][0])
                    z_y = np.append(z_y, XYZ[i][1])
                    z_z = np.append(z_z, XYZ[i][2])
                    V_z = np.append(V_z, VVV[i][0])
                else:
                    raise Exception()

            I_func['NS'] = NearestNDInterpolator((x_x, x_y, x_z), V_x)
            I_func['WE'] = NearestNDInterpolator((y_x, y_y, y_z), V_y)
            I_func['BF'] = NearestNDInterpolator((z_x, z_y, z_z), V_z)
        else:
            I_func = None
        I_func = cOmm.bcast(I_func, root=mAster_rank)
        func = (I_func['NS'], I_func['WE'], I_func['BF'])
        self.func._body_ = func
        self.___PRIVATE_discretize_the_flux_of_a_VectorField_of_ftype_standard___()
        self.func._body_ = None



if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\form\trace\_2_trace.py

    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy_periodic', c=0.0)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',7), ('Lobatto',8), ('Lobatto',9)])
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return t + np.cos(2*np.pi*x) * np.cos(4*np.pi*y) * np.cos(6*np.pi*z)
    S = FC('scalar', p)

    t2 = FC('2-t')

    # t2.TW.func.DO.set_func_body_as(S)
    # t2.TW.current_time = 0
    # t2.TW.DO.push_all_to_instant()

    # t2.discretize()
    #
    # xi = eta = sigma = np.linspace(-1,1,20)
    #
    # t2.visualize.matplot()

    # print(mesh.elements.map)
    # print(mesh.trace.elements.GLOBAL_num)
    # print(rAnk, t2.numbering.GLOBAL_boundary_dofs_ravel)
    # print(t2.coboundary.trace_matrix)

    # print(rAnk, xyz.keys())


    M = t2.___PRIVATE_generate_TEW_mass_matrices___()

    # print(M.keys())