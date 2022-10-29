# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/19 11:48 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np
from objects.CSCG._2d.mesh.domain.regions.region.interpolations.transfinite.mapping import TransfiniteMapping


class miUsGrid_TriangularMesh_Element_CT(FrozenOnly):
    """We map the LEFT edge of the reference element into the singular vertex.

    We use the transfinite mapping for the transformation.

    """
    def __init__(self, element):
        """"""
        self._element_ = element
        self._x0_, self._y0_ = self._element_.coordinates[0]
        self._x1_, self._y1_ = self._element_.coordinates[1]
        self._x2_, self._y2_ = self._element_.coordinates[2]
        self._x10 = self._x1_ - self._x0_
        self._y10 = self._y1_ - self._y0_
        self._x20 = self._x2_ - self._x0_
        self._y20 = self._y2_ - self._y0_
        self._x12 = self._x1_ - self._x2_
        self._y12 = self._y1_ - self._y2_
        self._TFM_ = TransfiniteMapping(
                (self._Pr_L_XY  , self._Pr_D_XY  , self._Pr_R_XY  , self._Pr_U_XY  ),
                (self._Pr_L_XoYo, self._Pr_D_XoYo, self._Pr_R_XoYo, self._Pr_U_XoYo)
        )
        self._freeze_self_()

    def _Pr_D_XY(self, o):
        return self._x0_ + o * self._x10, \
               self._y0_ + o * self._y10
    def _Pr_R_XY(self, o):
        return self._x2_ + o * self._x12, \
               self._y2_ + o * self._y12
    def _Pr_U_XY(self, o):
        return self._x0_ + o * self._x20, \
               self._y0_ + o * self._y20
    def _Pr_L_XY(self, o):
        return self._x0_ * np.ones_like(o), \
               self._y0_ * np.ones_like(o)

    def _Pr_D_XoYo(self, o):
        return self._x10 * np.ones_like(o), \
               self._y10 * np.ones_like(o)
    def _Pr_R_XoYo(self, o):
        return self._x12 * np.ones_like(o), \
               self._y12 * np.ones_like(o)
    def _Pr_U_XoYo(self, o):
        return self._x20 * np.ones_like(o), \
               self._y20 * np.ones_like(o)
    @staticmethod
    def _Pr_L_XoYo(o):
        ZO = np.zeros_like(o)
        return ZO, ZO



    def mapping(self, r, s):
        """ r, s be in [-1, 1]. """
        r = (r + 1) / 2
        s = (s + 1) / 2
        return self._TFM_.mapping(r, s)

    def Jacobian_matrix(self, r, s):
        """ r, s be in [-1, 1]. """
        r = (r + 1) / 2
        s = (s + 1) / 2
        return ((0.5 * self._TFM_.dx_dr(r, s), 0.5 * self._TFM_.dx_ds(r, s)),
                (0.5 * self._TFM_.dy_dr(r, s), 0.5 * self._TFM_.dy_ds(r, s)))

    def Jacobian(self, *evaluationPoints, itmD=None):
        """Determinant of the Jacobian matrix."""
        if itmD is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        else:
            J = itmD
        return J[0][0]*J[1][1] - J[0][1]*J[1][0]

    def metric(self, *evaluationPoints, itmD=None):
        """
        The metric ``g:= det(G):=(det(J))**2``. Since our Jacobian and inverse of Jacobian are both square,
        we know that the metric ``g`` is equal to square of ``det(J)``. ``g = (det(J))**2`` is due to the
        fact that the Jacobian matrix is square. The definition of ``g`` usually is given
        as ``g:= det(G)`` where ``G`` is the metric matrix, or metric tensor.
        """
        if itmD is None:
            detJ = self.Jacobian(*evaluationPoints)
        else:
            detJ = itmD
        return detJ ** 2

    def inverse_Jacobian_matrix(self, *evaluationPoints, itmD=None):
        """The inverse Jacobian matrix. """
        if itmD is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        else:
            J = itmD
        Jacobian = J[0][0]*J[1][1] - J[0][1]*J[1][0]
        reciprocalJacobian = 1 / Jacobian
        del Jacobian
        iJ00 = + reciprocalJacobian * J[1][1]
        iJ01 = - reciprocalJacobian * J[0][1]
        iJ10 = - reciprocalJacobian * J[1][0]
        iJ11 = + reciprocalJacobian * J[0][0]
        return [[iJ00, iJ01],
                [iJ10, iJ11]]

    def inverse_Jacobian(self, *evaluationPoints, itmD=None):
        """Determinant of the inverse Jacobian matrix. """
        if itmD is None:
            iJ = self.inverse_Jacobian_matrix(*evaluationPoints)
        else:
            iJ = itmD
        return iJ[0][0]*iJ[1][1] - iJ[0][1]*iJ[1][0]

    def metric_matrix(self, *evaluationPoints, itmD=None):
        """
        Also called metric tensor. Let J be the Jacobian matrix. The ``metricMatrix`` is
        denoted by G, G := J^T.dot(J). And the metric is ``g := (det(J))**2 or g := det(G).``
        Which means for a square Jacobian matrix, the metric turns out to be the square of the
        determinant of the Jacobian matrix.

        The entries of G is normally denoted as g_{i,j}.
        """
        if itmD is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        else:
            J = itmD
        G = [[None for _ in range(2)] for __ in range(2)]
        for i in range(2):
            for j in range(i, 2):
                # noinspection PyTypeChecker
                G[i][j] = J[0][i] * J[0][j]
                for l in range(1, 2):
                    G[i][j] += J[l][i] * J[l][j]
                if i != j:
                    G[j][i] = G[i][j]
        return G

    def inverse_metric_matrix(self, *evaluationPoints, itmD=None):
        """
        The ``inverseMetricMatrix`` is the metric matrix of the inverse Jacobian matrix
        or the metric of the inverse mapping. It is usually denoted as G^{-1}.

        The entries of G^{-1} is normally denoted as g^{i,j}.
        """
        if itmD is None:
            iJ = self.inverse_Jacobian_matrix(*evaluationPoints)
        else:
            iJ = itmD
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




if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/elements/element/coordinate_transformation/main.py
    pass
