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

from root.config.main import np
from screws.quadrature import Quadrature
from objects.CSCG._3d.forms.trace.base.main import _3dCSCG_Standard_Trace
from scipy.sparse import csr_matrix, bmat
from objects.CSCG._3d.forms.trace._1tr.discretize.main import _3dCSCG_1Trace_Discretize
from objects.CSCG._3d.forms.trace._1tr.visualize import _3dCSCG_1Trace_Visualize


class _3dCSCG_1Trace(_3dCSCG_Standard_Trace, ABC):
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
        self.___PRIVATE_reset_cache___()
        self._discretize_ = _3dCSCG_1Trace_Discretize(self)
        self._visualize_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard', 'trace-element-wise'), \
                f"3dCSCG 1-trace FUNC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."

        else:
            raise NotImplementedError(
                f"1-trace form cannot accommodate {func_body}.")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('trace-element-wise', ), \
                f"3dCSCG 1-trace BC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise NotImplementedError(
                f"3d CSCG 1-trace form BC cannot accommodate {func_body}.")

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_1Trace_Visualize(self)
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
                        v[key] = [v_x.reshape((jj, kk), order='F'),
                                  v_y.reshape((jj, kk), order='F'),
                                  v_z.reshape((jj, kk), order='F')]
                    elif side in 'WE':
                        xyz[key] = [xyz_i[m].reshape(ii, kk, order='F') for m in range(3)]
                        v[key] = [v_x.reshape((ii, kk), order='F'),
                                  v_y.reshape((ii, kk), order='F'),
                                  v_z.reshape((ii, kk), order='F')]
                    elif side in 'BF':
                        xyz[key] = [xyz_i[m].reshape(ii, jj, order='F') for m in range(3)]
                        v[key] = [v_x.reshape((ii, jj), order='F'),
                                  v_y.reshape((ii, jj), order='F'),
                                  v_z.reshape((ii, jj), order='F')]
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

                M00 = csr_matrix(M00)
                M11 = csr_matrix(M11)
                M01 = csr_matrix(M01)
                M10 = M01.T

                M = bmat([(M00, M01),
                          (M10, M11)], format='csr')

                if isinstance(mark, str): local_cache[mark] = M

                MD[i] = M

        return MD












if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\trace\_1_trace\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',4), ('Lobatto',4), ('Lobatto',4)])
    FC = FormCaller(mesh, space)

    def u(t, x, y, z): return t + np.sin(2*np.pi*x) * np.cos(np.pi*y) * np.cos(2*np.pi*z)
    def v(t, x, y, z): return t + np.cos(2*np.pi*x) * np.sin(2*np.pi*y) * np.cos(np.pi*z)
    def w(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.sin(2*np.pi*z)
    # def u(t, x, y, z): return 0.3 + 0*x
    # def v(t, x, y, z): return 0.75 + 0*x
    # def w(t, x, y, z): return 0.25 + 0*x
    V = FC('vector', (u,v,w))

    t1 = FC('1-t')

    t1.TW.func.do.set_func_body_as(V)
    t1.TW.current_time = 0
    t1.TW.do.push_all_to_instant()

    t1.discretize()

    t1.visualize()
    #
    # xi, et, sg = np.linspace(-1,1,9), np.linspace(-1,1,10), np.linspace(-1,1,11)
    #
    # xyz, v = t1.reconstruct(xi, et, sg)
    # print(v.keys())

    # t1.visualize()

    # tM1 = t1.matrices.mass





