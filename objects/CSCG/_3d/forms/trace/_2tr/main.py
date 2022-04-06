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
from root.save import read
from scipy.interpolate import NearestNDInterpolator
from screws.quadrature import Quadrature
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
        self.___PRIVATE_reset_cache___()
        self._discretize_ = _3dCSCG_2Trace_Discretize(self)
        self._visualize_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
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

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
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

    def ___PRIVATE_do_resemble___(self, obj_or_filename):
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
        self.discretize._standard_vector_()
        self.func._body_ = None






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