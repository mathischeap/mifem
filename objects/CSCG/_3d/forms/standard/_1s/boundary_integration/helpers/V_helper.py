# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/26/2022 10:38 AM
"""
import numpy as np

from components.freeze.main import FrozenOnly
from components.quadrature import Quadrature

from scipy.sparse import csc_matrix


class S1F_BI_V_Helper(FrozenOnly):
    """"""

    def __init__(self, s1f, V, quad_degree):
        """

        Parameters
        ----------
        s1f : standard-1-form
            The self 1-form. No need to have a cochain.
        V :
            A 3dCSCG vector.
        quad_degree
        """
        self._s1f_ = s1f
        self._V_ = V

        if quad_degree is None:
            quad_degree = s1f.dqp

        QUAD, WEIGHTs = Quadrature(quad_degree, category='Gauss').quad

        Qx, Qy, Qz = QUAD
        Wx, Wy, Wz = WEIGHTs

        Qn = [[-1, ], Qy, Qz]
        Qs = [[1, ], Qy, Qz]
        Qw = [Qx, [-1, ], Qz]
        Qe = [Qx, [1, ], Qz]
        Qb = [Qx, Qy, [-1, ]]
        Qf = [Qx, Qy, [1, ]]

        self.Wn = self.Ws = np.kron(Wz, Wy)
        self.Ww = self.We = np.kron(Wz, Wx)
        self.Wb = self.Wf = np.kron(Wy, Wx)

        self.GMn = s1f.do.make_reconstruction_matrix_on_grid(*Qn, element_range='mesh boundary')
        self.GMs = s1f.do.make_reconstruction_matrix_on_grid(*Qs, element_range='mesh boundary')
        self.GMw = s1f.do.make_reconstruction_matrix_on_grid(*Qw, element_range='mesh boundary')
        self.GMe = s1f.do.make_reconstruction_matrix_on_grid(*Qe, element_range='mesh boundary')
        self.GMb = s1f.do.make_reconstruction_matrix_on_grid(*Qb, element_range='mesh boundary')
        self.GMf = s1f.do.make_reconstruction_matrix_on_grid(*Qf, element_range='mesh boundary')

        self.Qn = np.meshgrid(*Qn, indexing='ij')
        self.Qs = np.meshgrid(*Qs, indexing='ij')
        self.Qw = np.meshgrid(*Qw, indexing='ij')
        self.Qe = np.meshgrid(*Qe, indexing='ij')
        self.Qb = np.meshgrid(*Qb, indexing='ij')
        self.Qf = np.meshgrid(*Qf, indexing='ij')

        self.mesh = s1f.mesh

        self.ZeroVec = csc_matrix(np.zeros((s1f.num.basis, 1)))

        self._freeze_self_()

    def __call__(self, i):
        """for ith mesh element."""
        element = self.mesh.elements[i]

        if element.whether.internal:
            return self.ZeroVec

        else:
            positions = element.position

            V = list()

            assert self._V_.current_time is not None, f"Set the current_time for the vector first!"
            t = self._V_.current_time

            for j, bn in enumerate(positions):

                if isinstance(bn, str):  # on mesh boundary
                    side = 'nswebf'[j]

                    vf0, vf1, vf2 = self._V_.func[bn]

                    vf0 = vf0(t, *getattr(self, 'Q'+side)).ravel('F')
                    vf1 = vf1(t, *getattr(self, 'Q'+side)).ravel('F')
                    vf2 = vf2(t, *getattr(self, 'Q'+side)).ravel('F')

                    GM0, GM1, GM2 = getattr(self, 'GM'+side)[i]

                    QW = getattr(self, 'W'+side)

                    Vs = \
                        np.einsum('w, wb -> wb', vf0, GM0, optimize='optimal') \
                        + np.einsum('w, wb -> wb', vf1, GM1, optimize='optimal') \
                        + np.einsum('w, wb -> wb', vf2, GM2, optimize='optimal')

                    TE = self.mesh.trace.elements.map[i][j]
                    TEct = self.mesh.trace.elements[TE].coordinate_transformation

                    C_Jacobian = TEct.constant.Jacobian

                    if C_Jacobian is not None:  # we have a constant Jacobian
                        Vs = csc_matrix(
                            np.einsum(
                                'wb, w -> b',
                                Vs,
                                QW * C_Jacobian,
                                optimize='optimal')[:, np.newaxis])

                    else:
                        raise NotImplementedError()

                    V.append(Vs)

                else:
                    pass

            V = np.sum(V, axis=0)
            return V
