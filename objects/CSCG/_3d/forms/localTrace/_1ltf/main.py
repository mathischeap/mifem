# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/26/2022 2:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.CSCG._3d.forms.localTrace.base.main import _3dCSCG_LocalTrace
from objects.CSCG._3d.forms.localTrace._1ltf.discretize.main import _3dCSCG_1LocalTrace_Discretize
from objects.CSCG._3d.forms.localTrace._1ltf.reconstruct import _3dCSCG_1LocalTrace_Reconstruct
from objects.CSCG._3d.forms.localTrace._1ltf.visualize.main import _3dCSCG_1LocalTrace_Visualize
from components.quadrature import Quadrature
import numpy as np
from scipy.sparse import csr_matrix, bmat

class _3dCSCG_1LocalTrace(_3dCSCG_LocalTrace):
    """"""

    def __init__(self, mesh, space, hybrid=True, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-1-local-trace-form'):
        super().__init__(mesh, space, hybrid, orientation, numbering_parameters, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_localtrace_1form')
        self._discretize_ = _3dCSCG_1LocalTrace_Discretize(self)
        self._reconstruct_ = _3dCSCG_1LocalTrace_Reconstruct(self)
        self._visualize_ = _3dCSCG_1LocalTrace_Visualize(self)
        self._freeze_self_()


    def ___Pr_check_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard', 'trace-element-wise'), \
                f"3dCSCG 1-ltf FUNC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."

        else:
            raise NotImplementedError(
                f"3d CSCG 1-ltf cannot accommodate {func_body}.")

    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('trace-element-wise', ), \
                f"3dCSCG 1-ltf BC cannot accommodate _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise NotImplementedError(
                f"3d CSCG 1-ltf BC cannot accommodate {func_body}.")

    def ___PrLT_mass_matrices___(self):
        """This is a brutal force version for the mesh-element-wise mass matrix."""
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

        R_sides = self.do.make_reconstruction_matrix_on_grid(*quad_nodes)[1]

        MD = dict()
        cacheDict = dict()

        Tmap = self.mesh.trace.elements.map
        for i in self.mesh.elements:

            element = self.mesh.elements[i]

            mark = element.type_wrt_metric.mark

            if isinstance(mark, str) and mark in cacheDict:
                MD[i] = cacheDict[mark]
            else:
                tes = Tmap[i]
                RSi = R_sides[i]
                M = list()
                for j, side in zip(tes, 'NSWEBF'):
                    te = self.mesh.trace.elements[j]
                    R = RSi[side] # array, if sparse, do toarray()

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

                    M.append(
                        csr_matrix(
                            np.einsum(
                                'vi, vj, v, v -> ij',
                                R, R, np.sqrt(g), QW,
                                optimize='optimal',
                            )
                        )
                    )

                if self.whether.hybrid:
                    M = bmat(
                        (
                            [M[0], None, None, None, None, None],
                            [None, M[1], None, None, None, None],
                            [None, None, M[2], None, None, None],
                            [None, None, None, M[3], None, None],
                            [None, None, None, None, M[4], None],
                            [None, None, None, None, None, M[5]]
                        ),
                        format='csr'
                    )
                else:
                    raise NotImplementedError()

                if isinstance(mark, str):
                    cacheDict[mark] = M
                else:
                    pass

                MD[i] = M

        return MD




if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/localTrace/_1ltf/main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller  # , ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([2, 2, 2])
    space = SpaceInvoker('polynomials')([3, 3, 3])
    FC = FormCaller(mesh, space)

    lt1 = FC('1-lt')
    print(lt1.numbering.gathering)