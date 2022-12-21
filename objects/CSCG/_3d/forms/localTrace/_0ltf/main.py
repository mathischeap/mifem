# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/26/2022 2:56 PM
"""
import sys
from abc import ABC

if './' not in sys.path:
    sys.path.append('./')
from objects.CSCG._3d.forms.localTrace.base.main import _3dCSCG_LocalTrace
from objects.CSCG._3d.forms.localTrace._0ltf.discretize.main import _3dCSCG_0LocalTrace_Discretize
from objects.CSCG._3d.forms.localTrace._0ltf.reconstruct import _3dCSCG_0LocalTrace_Reconstruct
from objects.CSCG._3d.forms.localTrace._0ltf.visualize.main import _3dCSCG_0LocalTrace_Visualize
from objects.CSCG._3d.forms.localTrace._0ltf.boundary_integration.main import _3dCSCG_0LTF_BoundaryIntegration
from components.quadrature import Quadrature
import numpy as np
from scipy.sparse import csr_matrix, bmat
from components.assemblers import MatrixAssembler


class _3dCSCG_0LocalTrace(_3dCSCG_LocalTrace, ABC):
    """"""

    def __init__(
        self, mesh, space, hybrid=True, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-0-local-trace-form',
    ):
        super().__init__(mesh, space, hybrid, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_localtrace_0form')
        self._discretize_ = _3dCSCG_0LocalTrace_Discretize(self)
        self._reconstruct_ = _3dCSCG_0LocalTrace_Reconstruct(self)
        self._visualize_ = _3dCSCG_0LocalTrace_Visualize(self)
        self._bi_ = None
        self._freeze_self_()

    def ___Pr_check_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard', 'boundary-wise'), \
                f"3dCSCG 0ltf FUNC does not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0ltf FUNC does not accept func {func_body.__class__}")

    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard', 'boundary-wise', 'trace-element-wise'), \
                f"3dCSCG 0ltf BC does not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0ltf BC does not accept func {func_body.__class__}")

    def ___PrLT_mass_matrices_brutal_force___(self):
        """Should return the same (mesh-element-wise) mass matrices as ___PrLT_mass_matrices___."""

        if self.mesh.whether.orthogonal:
            p = self.dqp
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
        else:
            p = [self.dqp[i] + 1 for i in range(self.ndim)]
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad

        x, y = np.meshgrid(quad_nodes[1], quad_nodes[2], indexing='ij')
        xNS = x.ravel('F')
        yNS = y.ravel('F')
        x, y = np.meshgrid(quad_nodes[0], quad_nodes[2], indexing='ij')
        xWE = x.ravel('F')
        yWE = y.ravel('F')
        x, y = np.meshgrid(quad_nodes[0], quad_nodes[1], indexing='ij')
        xBF = x.ravel('F')
        yBF = y.ravel('F')

        qw = dict()
        qw['NS'] = np.kron(quad_weights[2], quad_weights[1])
        qw['WE'] = np.kron(quad_weights[2], quad_weights[0])
        qw['BF'] = np.kron(quad_weights[1], quad_weights[0])

        R_sides = self.do.make_reconstruction_matrix_on_grid(*quad_nodes)

        MD = dict()
        cacheDict = dict()

        Tmap = self.mesh.trace.elements.map
        if self.whether.hybrid:
            assembler = None
        else:
            local_gathering = self.numbering.local_gathering
            assembler = MatrixAssembler(local_gathering, local_gathering)

        for i in self.mesh.elements:

            element = self.mesh.elements[i]

            mark = element.type_wrt_metric.mark

            if isinstance(mark, str) and mark in cacheDict:
                MD[i] = cacheDict[mark]
            else:
                tes = Tmap[i]
                RSi = R_sides[i]
                M = dict()
                for j, side in zip(tes, 'NSWEBF'):
                    te = self.mesh.trace.elements[j]
                    R = RSi[side]  # array, if sparse, do toarray()

                    if side in 'NS':
                        g = te.coordinate_transformation.metric(xNS, yNS)
                        QW = qw['NS']
                    elif side in 'WE':
                        g = te.coordinate_transformation.metric(xWE, yWE)
                        QW = qw['WE']
                    elif side in 'BF':
                        g = te.coordinate_transformation.metric(xBF, yBF)
                        QW = qw['BF']
                    else:
                        raise Exception()

                    M[side] = np.einsum(
                        'vi, vj, v, v -> ij',
                        R, R, np.sqrt(g), QW,
                        optimize='optimal',
                    )

                if self.whether.hybrid:
                    M = bmat(
                        (
                            [csr_matrix(M['N']), None, None, None, None, None],
                            [None, csr_matrix(M['S']), None, None, None, None],
                            [None, None, csr_matrix(M['W']), None, None, None],
                            [None, None, None, csr_matrix(M['E']), None, None],
                            [None, None, None, None, csr_matrix(M['B']), None],
                            [None, None, None, None, None, csr_matrix(M['F'])]
                        ),
                        format='csr'
                    )

                else:
                    M = assembler(M, 'add', format='csr')

                if isinstance(mark, str):
                    cacheDict[mark] = M
                else:
                    pass

                MD[i] = M

        return MD

    def ___PrLT_mass_matrices___(self):
        """Make the mesh-element-wise mass matrix."""
        if self.mesh.whether.orthogonal:
            p = self.dqp
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
        else:
            p = [self.dqp[i] + 1 for i in range(self.ndim)]
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad

        qw = dict()
        qw['NS'] = np.kron(quad_weights[2], quad_weights[1])
        qw['WE'] = np.kron(quad_weights[2], quad_weights[0])
        qw['BF'] = np.kron(quad_weights[1], quad_weights[0])

        xietasigma, pb = self.do.evaluate_basis_at_meshgrid(*quad_nodes)

        local_cache = dict()

        MD = dict()

        Tmap = self.mesh.trace.elements.map
        if self.whether.hybrid:
            assembler = None
        else:
            local_gathering = self.numbering.local_gathering
            assembler = MatrixAssembler(local_gathering, local_gathering)

        for i in self.mesh.elements:

            mark = self.mesh.elements[i].type_wrt_metric.mark

            if isinstance(mark, str) and mark in local_cache:
                MD[i] = local_cache[mark]
            else:
                TES = Tmap[i]
                M = dict()
                for j, side in zip(TES, 'NSWEBF'):
                    te = self.mesh.trace.elements[j]
                    g = te.coordinate_transformation.metric(*xietasigma[side])
                    b = pb[side][0]

                    if side in 'NS':
                        M_ = np.einsum(
                            'im, jm, m -> ij', b, b, np.sqrt(g) * qw['NS'], optimize='greedy'
                        )

                    elif side in 'WE':
                        M_ = np.einsum(
                            'im, jm, m -> ij', b, b, np.sqrt(g) * qw['WE'], optimize='greedy'
                        )

                    elif side in 'BF':
                        M_ = np.einsum(
                            'im, jm, m -> ij', b, b, np.sqrt(g) * qw['BF'], optimize='greedy'
                        )

                    else:
                        raise Exception()

                    M[side] = M_

                if self.whether.hybrid:
                    M = bmat(
                        (
                            [csr_matrix(M['N']), None, None, None, None, None],
                            [None, csr_matrix(M['S']), None, None, None, None],
                            [None, None, csr_matrix(M['W']), None, None, None],
                            [None, None, None, csr_matrix(M['E']), None, None],
                            [None, None, None, None, csr_matrix(M['B']), None],
                            [None, None, None, None, None, csr_matrix(M['F'])]
                        ),
                        format='csr'
                    )
                else:
                    M = assembler(M, 'add', format='csr')

                if isinstance(mark, str):
                    local_cache[mark] = M
                else:
                    pass

                MD[i] = M

        return MD

    @property
    def boundary_integration(self):
        if self._bi_ is None:
            self._bi_ = _3dCSCG_0LTF_BoundaryIntegration(self)
        return self._bi_


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/localTrace/_0ltf/main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([3, 3, 3])
    space = SpaceInvoker('polynomials')([1, 1, 1])
    FC = FormCaller(mesh, space)

    ltf0 = FC('0-lt', hybrid=False)

    def p(t, x, y, z): return np.sin(2*np.pi*x) + t + 0 * y * z
    scalar = FC('scalar', {'North': p, 'South': p, 'West': p, 'East': p, 'Back': p, 'Front': p})


    ltf0.BC.CF = scalar
    ltf0.BC.boundaries = 'all'
    scalar.current_time = 0

    BI = ltf0.BC.interpret.boundary_integration

    for i in BI:
        _ = BI[i]
        # print(i, BI[i])
