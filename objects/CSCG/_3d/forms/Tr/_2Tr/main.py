

import sys
if './' not in sys.path: sys.path.append('./')

import numpy as np
from screws.quadrature import Quadrature

from objects.CSCG._3d.forms.Tr.base.main import _3dCSCG_Standard_Tr
from objects.CSCG._3d.forms.Tr._2Tr.discretize.main import _3dCSCG_2Tr_Discretize
from objects.CSCG._3d.forms.Tr._2Tr.visualize import _3dCSCG_Tr_Visualize
from objects.CSCG._3d.forms.Tr._2Tr.reconstruct import _3dCSCG_2Tr_Reconstruct




class _3dCSCG_2Tr(_3dCSCG_Standard_Tr):
    """
    Tr 2-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-2-Tr-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_Tr_2form')
        self.___PRIVATE_reset_cache___()
        self._discretize_ = None
        self._reconstruct_ = None
        self._visualize_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super(_3dCSCG_2Tr, self).___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 2-Tr FUNC cannot accommodate _3dCSCG_ScalarField of ftype {func_body.ftype}."
        elif func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 2-Tr FUNC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise NotImplementedError(
                f"3d CSCG 2-Tr form FUNC cannot accommodate {func_body}.")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard', 'boundary-wise'), \
                f"3dCSCG 2-Tr BC cannot accommodate _3dCSCG_ScalarField of ftype {func_body.ftype}."
        elif func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard', 'boundary-wise'), \
                f"3dCSCG 2-Tr BC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise NotImplementedError(
                f"3d CSCG 2-Tr form BC cannot accommodate {func_body}.")

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_Tr_Visualize(self)
        return self._visualize_

    @property
    def discretize(self):
        if self._discretize_ is None:
            self._discretize_ = _3dCSCG_2Tr_Discretize(self)
        return self._discretize_

    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = _3dCSCG_2Tr_Reconstruct(self)
        return self._reconstruct_


    def ___PRIVATE_generate_TEW_mass_matrices___(self):
        """Generate the trace-element-wise mass matrices."""

        p = [self.dqp[i] for i in range(self.ndim)]
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
                    M = np.einsum('im, jm, m -> ij',
                                  b, b, np.reciprocal(np.sqrt(g)) * qw['NS'], optimize='greedy')
                elif side in 'WE':
                    M = np.einsum('im, jm, m -> ij',
                                  b, b, np.reciprocal(np.sqrt(g)) * qw['WE'], optimize='greedy')
                elif side in 'BF':
                    M = np.einsum('im, jm, m -> ij',
                                  b, b, np.reciprocal(np.sqrt(g)) * qw['BF'], optimize='greedy')
                else:
                    raise Exception()

                if isinstance(mark, str): local_cache[mark] = M

                MD[i] = M

        return MD






if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\Tr\_2Tr\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, \
        FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([1, 1, 1])
    space = SpaceInvoker('polynomials')([('Lobatto', 5), ('Lobatto', 5), ('Lobatto', 5)])
    FC = FormCaller(mesh, space)
    T2 = FC('2-Tr')
    def p(t, x, y, z): return t + np.cos(2*np.pi*x) * np.cos(2*np.pi*y) * np.cos(2*np.pi*z)
    S = FC('scalar', p)
    T2.TW.func.body = S
    T2.TW.current_time = 0
    T2.TW.do.push_all_to_instant()
    T2.discretize()

    T = T2.matrices.trace
    S = T2.matrices.selective
    M = T2.matrices.mass_TEW