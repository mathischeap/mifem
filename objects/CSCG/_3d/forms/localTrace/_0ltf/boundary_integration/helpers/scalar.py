# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/19/2022 1:41 PM
"""
from components.freeze.main import FrozenOnly
from scipy.sparse import csc_matrix
from components.quadrature import Quadrature
import numpy as np


class Scalar(FrozenOnly):
    """"""

    def __init__(self, ltf, scalar, quad_degree):
        """"""
        self._ltf_ = ltf
        self._scalar_ = scalar
        quad = Quadrature(quad_degree, category='Gauss')

        nodes, weights = quad.quad
        Wx, Wy, Wz = weights
        self.Wn = self.Ws = np.kron(Wz, Wy)
        self.Ww = self.We = np.kron(Wz, Wx)
        self.Wb = self.Wf = np.kron(Wy, Wx)

        self.rm = self._ltf_.reconstruct.___PrLT_make_reconstruction_matrix_on_grid___(
            *nodes
        )[1]

        self.ZeroVec = csc_matrix((ltf.num.basis, 1))
        self.mesh = ltf.mesh
        self._freeze_self_()

    def __call__(self, i):
        """"""
        element = self.mesh.elements[i]
        if element.whether.internal:
            return self.ZeroVec

        else:
            raise NotImplementedError()
