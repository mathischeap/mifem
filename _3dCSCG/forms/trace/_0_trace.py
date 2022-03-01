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
from _3dCSCG.forms.trace.base.main import _3dCSCG_Standard_Trace
from screws.quadrature import Quadrature


class _3dCSCG_0Trace(_3dCSCG_Standard_Trace, ABC):
    """
    Trace 0-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-0-trace-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_trace_0form')
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0-trace FUNC cannot accommodate _3dCSCG_ScalarField of ftype {func_body.ftype}."
        elif func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0-trace FUNC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise NotImplementedError(
                f"0-trace form cannot accommodate {func_body}.")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0-trace BC cannot accommodate _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise NotImplementedError(
                f"3d CSCG 0-trace form BC cannot accommodate {func_body}.")


    def discretize(self, update_cochain=True, target='func'):
        """
        Do the discretization.

        :param bool update_cochain: Whether we update the cochain if the trace form.
        :param str target:
        :return: The cochain corresponding to the particular discretization scheme.
        """
        if target == 'func':
            if self.TW.func.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(
                        target = 'func',
                        update_cochain=update_cochain)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target func) discretize _3dCSCG_ScalarField of ftype {self.func.ftype}.')

            elif self.TW.func.body.__class__.__name__ == '_3dCSCG_VectorField':
                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_the_flux_of_a_VectorField_of_standard_ftype___(
                        target='func',
                        update_cochain=update_cochain)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target func) discretize _3dCSCG_VectorField of ftype {self.func.ftype}.')

            else:
                raise NotImplementedError(f'3dCSCG 0-trace can not (target func) discretize {self.TW.func.body.__class__}.')


        elif target == 'BC':
            if self.TW.BC.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if self.BC.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(
                        target = 'BC',
                        update_cochain=False)
                else:
                    raise Exception(f'3dCSCG 0-trace can not (target BC) discretize _3dCSCG_ScalarField of ftype {self.BC.ftype}.')
            else:
                raise NotImplementedError(f'3dCSCG 0-trace can not (target BC) discretize {self.TW.BC.body.__class__}.')
        else:
            raise NotImplementedError()

    def ___PRIVATE_discretize_standard_ftype___(self, target='func', update_cochain=True):
        """We will discretize the standard scalar field to all trace elements.

        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.

        """
        if target in ('BC',): assert update_cochain is False, f"CANNOT update cochain when target is {target}"

        nodes = self.space.nodes
        nx, ny, nz = nodes
        nodes_NS = np.meshgrid(ny, nz, indexing='ij')
        nodes_WE = np.meshgrid(nx, nz, indexing='ij')
        nodes_BF = np.meshgrid(nx, ny, indexing='ij')

        if target == 'func':
            assert self.func.body is not None, f"No func.body!"
            _lf_ = self.func.body[0]
        elif target == 'BC':
            assert self.BC.body is not None, f"No BC.body!"
            _lf_ = self.BC.body[0]
        else:
            raise NotImplementedError(f"Not applicable for target={target}.")

        local_TEW = dict()
        for i in self.mesh.trace.elements:
            te = self.mesh.trace.elements[i]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side
            if ele_side in 'NS':
                x, y, z = te.coordinate_transformation.mapping(*nodes_NS, from_element=ele, side=ele_side)
            elif ele_side in 'WE':
                x, y, z = te.coordinate_transformation.mapping(*nodes_WE, from_element=ele, side=ele_side)
            elif ele_side in 'BF':
                x, y, z = te.coordinate_transformation.mapping(*nodes_BF, from_element=ele, side=ele_side)
            else:
                raise Exception()

            te_primal_local = _lf_(x, y, z).ravel('F')

            if not self.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[i] = te_primal_local

        if update_cochain: self.cochain.local_TEW = local_TEW

        return 'locally full local TEW cochain', local_TEW

    def ___PRIVATE_discretize_the_flux_of_a_VectorField_of_standard_ftype___(self, target='func', update_cochain=True):
        """We will discretize the standard vector field (norm flux) to all trace elements.

        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.
        """
        if target in ('BC',): assert update_cochain is False, f"CANNOT update cochain when target is {target}"

        nodes = self.space.nodes
        nx, ny, nz = nodes
        nodes_NS = np.meshgrid(ny, nz, indexing='ij')
        nodes_WE = np.meshgrid(nx, nz, indexing='ij')
        nodes_BF = np.meshgrid(nx, ny, indexing='ij')

        if target == 'func':
            assert self.func.body is not None, f"No func.body!"
            fx, fy, fz = self.func.body
        else:
            raise NotImplementedError(f"Not applicable for target={target}.")

        local_TEW = dict()
        for i in self.mesh.trace.elements:
            te = self.mesh.trace.elements[i]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side
            if ele_side in 'NS':
                x, y, z = te.coordinate_transformation.mapping(*nodes_NS, from_element=ele, side=ele_side)
                ouv = te.coordinate_transformation.unit_normal_vector(*nodes_NS)

            elif ele_side in 'WE':
                x, y, z = te.coordinate_transformation.mapping(*nodes_WE, from_element=ele, side=ele_side)
                ouv = te.coordinate_transformation.unit_normal_vector(*nodes_WE)
            elif ele_side in 'BF':
                x, y, z = te.coordinate_transformation.mapping(*nodes_BF, from_element=ele, side=ele_side)
                ouv = te.coordinate_transformation.unit_normal_vector(*nodes_BF)
            else:
                raise Exception()

            f0, f1, f2 = fx(x, y, z), fy(x, y, z), fz(x, y, z)

            f = f0 * ouv[0] + f1 * ouv[1] + f2 * ouv[2]

            te_primal_local = f.ravel('F')

            if not self.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[i] = te_primal_local

        if update_cochain: self.cochain.local_TEW = local_TEW

        return 'locally full local TEW cochain', local_TEW





    def reconstruct(self, xi, eta, sigma, ravel=False, i=None):
        """Do the reconstruction.

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
                indices = [i, ]
            else:
                indices = i

        xietasigma, pb = self.do.evaluate_basis_at_meshgrid(
            xi, eta,sigma)
        ii, jj, kk = np.size(xi), np.size(eta), np.size(sigma)
        xyz = dict()
        v = dict()
        for i in indices:
            if i in self.mesh.trace.elements:
                te = self.mesh.trace.elements[i]
                ele = te.CHARACTERISTIC_element
                side = te.CHARACTERISTIC_side
                xyz_i = te.coordinate_transformation.mapping(
                    *xietasigma[side], from_element=ele, side=side)
                prime_cochain = self.cochain.local_TEW[i]

                vi = np.einsum('i, ij -> j', prime_cochain,
                               pb[side][0], optimize='greedy')

                if ravel:
                    xyz[i] = xyz_i
                    v[i] = [vi, ]
                else:
                    if side in 'NS':
                        xyz[i] = [xyz_i[m].reshape(jj, kk, order='F') for m in range(3)]
                        v[i] = [vi.reshape((jj, kk), order='F'), ]
                    elif side in 'WE':
                        xyz[i] = [xyz_i[m].reshape(ii, kk, order='F') for m in range(3)]
                        v[i] = [vi.reshape((ii, kk), order='F'), ]
                    elif side in 'BF':
                        xyz[i] = [xyz_i[m].reshape(ii, jj, order='F') for m in range(3)]
                        v[i] = [vi.reshape((ii, jj), order='F'), ]
                    else:
                        raise Exception
        return xyz, v




    def ___PRIVATE_generate_TEW_mass_matrices___(self):
        """Generate the trace-element-wise mass matrices stored in a dict whose keys are trace-element numbers
        and values are the mass matrices in the corresponding trace-elements.
        """
        p = [self.dqp[i]+3 for i in range(self.ndim)] # +3 for safety, the mass matrices of standard forms use dqp
        quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad

        qw = dict()
        qw['NS'] = np.kron(quad_weights[2], quad_weights[1])
        qw['WE'] = np.kron(quad_weights[2], quad_weights[0])
        qw['BF'] = np.kron(quad_weights[1], quad_weights[0])

        xietasigma, pb = self.do.evaluate_basis_at_meshgrid(*quad_nodes)

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
                    M = np.einsum('im, jm, m -> ij', b, b, np.sqrt(g) * qw['NS'], optimize='greedy')
                elif side in 'WE':
                    M = np.einsum('im, jm, m -> ij', b, b, np.sqrt(g) * qw['WE'], optimize='greedy')
                elif side in 'BF':
                    M = np.einsum('im, jm, m -> ij', b, b, np.sqrt(g) * qw['BF'], optimize='greedy')
                else:
                    raise Exception()

                if isinstance(mark, str): local_cache[mark] = M

                MD[i] = M

        return MD






if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\form\trace\_0_trace.py

    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)

    from numpy import sin, pi, cos

    def flux_func(t, x, y, z):
        return cos(2*pi*x) * cos(2*pi*y) * cos(2*pi*z) + t

    flux = FC('scalar', flux_func)
    t0 = FC('0-t')

    # t0.TW.func.do.set_func_body_as(flux)
    # t0.TW.current_time = 1
    # t0.TW.do.push_all_to_instant()
    # t0.discretize()


    def u(t, x, y, z): return t + np.cos(2*np.pi*x) * np.cos(np.pi*y) * np.cos(2*np.pi*z)
    def v(t, x, y, z): return t + np.cos(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(np.pi*z)
    def w(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.cos(2*np.pi*z)
    V = FC('vector', (u,v,w))

    # t0.TW.func.do.set_func_body_as(V)
    # t0.TW.current_time = 1
    # t0.TW.do.push_all_to_instant()
    # t0.discretize()
    #
    # t0.visualize()

    M = t0.matrices.mass