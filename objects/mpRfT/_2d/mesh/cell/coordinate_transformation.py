# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from screws.freeze.base import FrozenOnly

class mpRfT2_Mesh_Cell_CT(FrozenOnly):
    """"""
    def __init__(self, cell):
        self._cell_ = cell
        self._me_ct_ = self._cell_.mesh.cscg.elements[self._cell_.indices[0]].coordinate_transformation
        self._freeze_self_()

    @property
    def origin_and_delta(self):
        """consider the level-0 cell (cscg mesh-element) as a [-1,1]^2 reference domain.
        """
        return self._cell_.mesh.do.find.origin_and_delta(self._cell_.indices)

    def ___PRIVATE_plot_data___(self, density=None, zoom=1):
        """

        Parameters
        ----------
        density
        zoom

        Returns
        -------
        data : ndarray
            of shape (4, 3, density) which stands for 4 lines, xyz coordinates and density.

        """
        if density is None:
            mark = self._cell_.type_wrt_metric.mark
            if isinstance(mark, str) and mark[:4] == 'Orth':
                density = 2
            else:
                density = 10
        else:
            pass

        assert 0 < zoom <= 1, f"zoom={zoom} is wrong!"

        O = np.ones(density) * zoom
        M = - np.ones(density) * zoom
        S = np.linspace(-1 * zoom, 1 * zoom, density)

        data = np.array([
        np.array(self.mapping(M, S)), # U
        np.array(self.mapping(O, S)), # D
        np.array(self.mapping(S, M)), # L
        np.array(self.mapping(S, O)), # R
        ])

        return data

    def mapping(self, xi, et):
        """

        Parameters
        ----------
        xi :
            Any dimension, in [-1,1]. xi, et same shape.
        et :
            Any dimension, in [-1,1]. xi, et same shape.

        Returns
        -------

        """
        if self._cell_.level == 0:
            return self._me_ct_.mapping(xi, et)

        o, d = self.origin_and_delta
        xi = o[0] + (xi + 1) * d / 2
        et = o[1] + (et + 1) * d / 2
        return self._me_ct_.mapping(xi, et)

    def Jacobian(self, xi, et):
        """

        Parameters
        ----------
        xi :
            Any dimension, in [-1,1]. xi, et same shape.
        et :
            Any dimension, in [-1,1]. xi, et same shape.

        Returns
        -------

        """
        if self._cell_.level == 0:
            return self._me_ct_.Jacobian(xi, et)

        o, d = self.origin_and_delta
        xi = o[0] + (xi + 1) * d / 2
        et = o[1] + (et + 1) * d / 2
        return self._me_ct_.Jacobian(xi, et) * d**2 / 4

    def inverse_Jacobian(self, xi, et, iJ=None):
        """Determinant of the inverse Jacobian matrix. """
        if iJ is None:
            iJ = self.inverse_Jacobian_matrix(xi, et)
        return iJ[0][0]*iJ[1][1] - iJ[0][1]*iJ[1][0]


    def Jacobian_matrix(self, xi, et):
        """

        Parameters
        ----------
        xi :
            Any dimension, in [-1,1]. xi, et same shape.
        et :
            Any dimension, in [-1,1]. xi, et same shape.

        Returns
        -------

        """
        if self._cell_.level == 0:
            return self._me_ct_.Jacobian_matrix(xi, et)

        o, d = self.origin_and_delta
        xi = o[0] + (xi + 1) * d / 2
        et = o[1] + (et + 1) * d / 2
        J00_01, J10_11 = self._me_ct_.Jacobian_matrix(xi, et)
        J00, J01 = J00_01
        J10, J11 = J10_11

        return ([J00 * d / 2, J01 * d / 2],
                [J10 * d / 2, J11 * d / 2])


    def inverse_Jacobian_matrix(self, xi, et, J=None):
        """The inverse Jacobian matrix. """
        if J is None:
            J = self.Jacobian_matrix(xi, et)
        Jacobian = J[0][0]*J[1][1] - J[0][1]*J[1][0]
        reciprocalJacobian = 1 / Jacobian
        del Jacobian
        iJ00 = + reciprocalJacobian * J[1][1]
        iJ01 = - reciprocalJacobian * J[0][1]
        iJ10 = - reciprocalJacobian * J[1][0]
        iJ11 = + reciprocalJacobian * J[0][0]
        return [[iJ00, iJ01],
                [iJ10, iJ11]]



    def inverse_metric_matrix(self, *evaluationPoints, iJ=None):
        """
        The ``inverseMetricMatrix`` is the metric matrix of the inverse Jacobian matrix
        or the metric of the inverse mapping. It is usually denoted as G^{-1}.

        The entries of G^{-1} is normally denoted as g^{i,j}.
        """
        if iJ is None:
            iJ = self.inverse_Jacobian_matrix(*evaluationPoints)
        iG = [[None for _ in range(2)] for __ in range(2)]
        for i in range(2):
            for j in range(i, 2):
                # noinspection PyTypeChecker
                iG[i][j] = iJ[i][0] * iJ[j][0]
                for l in range(1, 2):
                    iG[i][j] += iJ[i][l] * iJ[j][l]
                if i != j:
                    iG[j][i] = iG[i][j]
        return iG



if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/cell/coordinate_transformation.py
    pass