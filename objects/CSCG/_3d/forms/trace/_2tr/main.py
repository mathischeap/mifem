# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from components.quadrature import Quadrature
from objects.CSCG._3d.forms.trace.base.main import _3dCSCG_Standard_Trace
from objects.CSCG._3d.forms.trace._2tr.discretize.main import _3dCSCG_2Trace_Discretize
from objects.CSCG._3d.forms.trace._2tr.visualize import _3dCSCG_2Trace_Visualize


class _3dCSCG_2Trace(_3dCSCG_Standard_Trace):
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
        self._discretize_ = _3dCSCG_2Trace_Discretize(self)
        self._visualize_ = None
        self._freeze_self_()

    def ___Pr_check_CF___(self, func_body):
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

    def ___Pr_check_BC_CF___(self, func_body):
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

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_2Trace_Visualize(self)
        return self._visualize_

    @property
    def discretize(self):
        return self._discretize_

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
        xietasigma, pb = self.do.evaluate_basis_at_meshgrid(xi, eta, sigma)
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






if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\trace\_2_trace\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return t + np.cos(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(2*np.pi*z)
    S = FC('scalar', p)

    t2 = FC('2-t')
    t2.TW.func.body = S
    t2.TW.current_time = 0
    t2.TW.do.push_all_to_instant()
    t2.discretize()

    T = t2.matrices.trace
    S = t2.matrices.selective
    M = t2.matrices.mass

    t2.visualize()