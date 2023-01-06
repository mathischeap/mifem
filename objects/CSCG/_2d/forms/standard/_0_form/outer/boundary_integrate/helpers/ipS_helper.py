# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 7:34 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')

from components.freeze.base import FrozenOnly
from components.quadrature import Quadrature
import numpy as np
from scipy.sparse import csc_matrix


class ipS_Helper(FrozenOnly):
    """"""

    def __init__(self, sf, S, quad_degree):
        """

        Parameters
        ----------
        sf : standard-0-form
            The self 0-form. No need to have a cochain.
        S :
            A 2dCSCG scalar.
        quad_degree :


        """
        self._sf_ = sf
        self._S_ = S

        if quad_degree is None:
            quad_degree = sf.dqp

        QUAD, WEIGHTs = Quadrature(quad_degree, category='Gauss').quad

        Qud, Qlr = QUAD
        Wud, Wlr = WEIGHTs

        Qu = [[-1, ], Qlr]
        Qd = [[1, ], Qlr]
        Ql = [Qud, [-1, ]]
        Qr = [Qud, [1, ]]
        self.Wu = self.Wd = Wud
        self.Wl = self.Wr = Wlr

        self.GMu = sf.do.make_reconstruction_matrix_on_grid(*Qu, element_range='mesh boundary')
        self.GMd = sf.do.make_reconstruction_matrix_on_grid(*Qd, element_range='mesh boundary')
        self.GMl = sf.do.make_reconstruction_matrix_on_grid(*Ql, element_range='mesh boundary')
        self.GMr = sf.do.make_reconstruction_matrix_on_grid(*Qr, element_range='mesh boundary')

        self.Qu = np.meshgrid(*Qu, indexing='ij')
        self.Qd = np.meshgrid(*Qd, indexing='ij')
        self.Ql = np.meshgrid(*Ql, indexing='ij')
        self.Qr = np.meshgrid(*Qr, indexing='ij')

        self.mesh = sf.mesh
        self.ZeroVec = csc_matrix(np.zeros((sf.num.basis, 1)))
        self._freeze_self_()



    def __call__(self, i):
        """for ith mesh element."""
        element = self.mesh.elements[i]

        if element.whether.internal:
            return self.ZeroVec

        else:
            positions = element.position

            V = list()

            assert self._S_.current_time is not None, f"Set the current_time for the vector first!"

            for j, bn in enumerate(positions):

                if isinstance(bn, str):  # on mesh boundary
                    side = 'udlr'[j]

                    fv = self._S_.___DO_evaluate_func_at_time___()[bn][0]

                    fv = fv(*getattr(self, 'Q'+side)).ravel('F')
                    GM = getattr(self, 'GM'+side)[i]
                    QW = getattr(self, 'W' + side)

                    TE = self.mesh.trace.elements.map[i][j]
                    TEct = self.mesh.trace.elements[TE].coordinate_transformation
                    C_Jacobian = TEct.constant.Jacobian


                    if C_Jacobian is not None:
                        Vs = csc_matrix(
                                np.einsum(
                                        'w, wb, w -> b',
                                        fv, GM,  QW * C_Jacobian,
                                        optimize='optimal')[:, np.newaxis])

                    else:
                        raise NotImplementedError()
                    V.append(Vs)

                else:
                    pass

            return sum(V)


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/_0_form/outer/boundary_integrate/helpers/ipS_helper.py
    pass
