# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/26/2022 10:38 AM
"""
import sys

import numpy as np

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from screws.quadrature import Quadrature


class S1F_BI_TpHelper(FrozenOnly):
    """"""

    def __init__(self, s1f, VTP, quad_degree):
        """"""
        self._s1f_ = s1f
        self._VTP_ = VTP

        if quad_degree is None:
            quad_degree = s1f.dqp

        QUAD, WEIGHTs = Quadrature(quad_degree, category='Gauss').quad

        Qx, Qy, Qz = QUAD
        Wx, Wy, Wz = WEIGHTs

        self.Qn = [[-1,], Qy, Qz]
        self.Qs = [[ 1,], Qy, Qz]
        self.Qw = [Qx, [-1,], Qz]
        self.Qe = [Qx, [ 1,], Qz]
        self.Qb = [Qx, Qy, [-1,]]
        self.Qf = [Qx, Qy, [ 1,]]

        self.Wns = np.kron(Wz, Wy)
        self.Wwe = np.kron(Wz, Wx)
        self.Wbf = np.kron(Wy, Wx)

        self.GMn = s1f.do.make_reconstruction_matrix_on_grid(*self.Qn)
        self.GMs = s1f.do.make_reconstruction_matrix_on_grid(*self.Qs)
        self.GMw = s1f.do.make_reconstruction_matrix_on_grid(*self.Qw)
        self.GMe = s1f.do.make_reconstruction_matrix_on_grid(*self.Qe)
        self.GMb = s1f.do.make_reconstruction_matrix_on_grid(*self.Qb)
        self.GMf = s1f.do.make_reconstruction_matrix_on_grid(*self.Qf)

        self.Qn_mg = np.meshgrid(*self.Qn,  indexing='ij')
        self.Qs_mg = np.meshgrid(*self.Qs,  indexing='ij')
        self.Qw_mg = np.meshgrid(*self.Qw,  indexing='ij')
        self.Qe_mg = np.meshgrid(*self.Qe,  indexing='ij')
        self.Qb_mg = np.meshgrid(*self.Qb,  indexing='ij')
        self.Qf_mg = np.meshgrid(*self.Qf,  indexing='ij')

        self.mesh = s1f.mesh

        self._freeze_self_()

    def __call__(self, basic_unit):
        """"""
        element = self.mesh.elements[basic_unit]

        # vfN = self._VTP_











if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass