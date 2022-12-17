# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly
import numpy as np


class _3dCSCG_Mesh_Element_CT(FrozenOnly):
    """The Coordinate transformation class for the mesh element."""
    def __init__(self, element):
        self._element_ = element
        self._mesh_ = self._element_._elements_._mesh_
        self._region_ = self._mesh_.domain.regions[self._element_.in_region]
        self._origin_ = None
        self._delta_ = None
        self._freeze_self_()

    @property
    def origin(self):
        if self._origin_ is None:
            in_region, local_indices = self._mesh_.do.find.region_name_and_local_indices_of_element(
                self._element_.i
            )
            self._origin_, self._delta_ = \
                self._mesh_.do.find.reference_origin_and_size_of_element_of_given_local_indices(
                    in_region, local_indices
                )
        return self._origin_

    @property
    def delta(self):
        if self._delta_ is None:
            in_region, local_indices = self._mesh_.do.find.region_name_and_local_indices_of_element(
                self._element_.i
            )
            self._origin_, self._delta_ = \
                self._mesh_.do.find.reference_origin_and_size_of_element_of_given_local_indices(
                    in_region, local_indices
                )
        return self._delta_

    def mapping(self, *evaluationPoints):
        XYZ = self._region_.interpolation(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return XYZ

    def X(self, *evaluationPoints):
        X = self._region_.interpolation.mapping_X(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return X

    def Y(self, *evaluationPoints):
        Y = self._region_.interpolation.mapping_Y(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Y

    def Z(self, *evaluationPoints):
        Z = self._region_.interpolation.mapping_Z(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Z

    def Jacobian_matrix(self, *evaluationPoints):
        mark = self._element_.type_wrt_metric.mark
        if isinstance(mark, str) and mark[:4] == 'Orth':
            xyz_xietasigma = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            rst = [(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)]
            J00 = self._region_.interpolation.Jacobian_Xr(*rst)
            J11 = self._region_.interpolation.Jacobian_Ys(*rst)
            J22 = self._region_.interpolation.Jacobian_Zt(*rst)
            xyz_xietasigma[0][0] = J00 * (self.delta[0] / 2)
            xyz_xietasigma[1][1] = J11 * (self.delta[1] / 2)
            xyz_xietasigma[2][2] = J22 * (self.delta[2] / 2)
        else:
            xyz_xietasigma = [[np.zeros(np.shape(evaluationPoints[j])) for _ in range(3)] for j in range(3)]
            xyz_rst = self._region_.interpolation.Jacobian_matrix(
                *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
            for j in range(3):
                for L in range(3):
                    xyz_xietasigma[j][L] = xyz_rst[j][L] * (self.delta[L] / 2)
        return xyz_xietasigma


    def J00(self, *evaluationPoints):
        Xr = self._region_.interpolation.Jacobian_Xr(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Xr * self.delta[0] / 2

    def J01(self, *evaluationPoints):
        Xs = self._region_.interpolation.Jacobian_Xs(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Xs * self.delta[1] / 2

    def J02(self, *evaluationPoints):
        Xt = self._region_.interpolation.Jacobian_Xt(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Xt * self.delta[2] / 2

    def J10(self, *evaluationPoints):
        Yr = self._region_.interpolation.Jacobian_Yr(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Yr * self.delta[0] / 2

    def J11(self, *evaluationPoints):
        Ys = self._region_.interpolation.Jacobian_Ys(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Ys * self.delta[1] / 2

    def J12(self, *evaluationPoints):
        Yt = self._region_.interpolation.Jacobian_Yt(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Yt * self.delta[2] / 2

    def J20(self, *evaluationPoints):
        Zr = self._region_.interpolation.Jacobian_Zr(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Zr * self.delta[0] / 2

    def J21(self, *evaluationPoints):
        Zs = self._region_.interpolation.Jacobian_Zs(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Zs * self.delta[1] / 2

    def J22(self, *evaluationPoints):
        Zt = self._region_.interpolation.Jacobian_Zt(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(3)])
        return Zt * self.delta[2] / 2

    def J0_(self, *evaluationPoints):
        return self.J00(*evaluationPoints), self.J01(*evaluationPoints), self.J02(*evaluationPoints)

    def J1_(self, *evaluationPoints):
        return self.J10(*evaluationPoints), self.J11(*evaluationPoints), self.J12(*evaluationPoints)

    def J2_(self, *evaluationPoints):
        return self.J20(*evaluationPoints), self.J21(*evaluationPoints), self.J22(*evaluationPoints)

    def Jacobian(self, *evaluationPoints, J=None):
        """Determinant of the Jacobian matrix."""
        if J is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        return + J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][0] \
            + J[0][2]*J[1][0]*J[2][1] - J[0][0]*J[1][2]*J[2][1] \
            - J[0][1]*J[1][0]*J[2][2] - J[0][2]*J[1][1]*J[2][0]

    def metric(self, *evaluationPoints, detJ=None):
        """
        The metric ``g:= det(G):=(det(J))**2``. Since our Jacobian and inverse of Jacobian are both square,
        we know that the metric ``g`` is equal to square of ``det(J)``. ``g = (det(J))**2`` is due to the
        fact that the Jacobian matrix is square. The definition of ``g`` usually is given
        as ``g:= det(G)`` where ``G`` is the metric matrix, or metric tensor.
        """
        if detJ is None:
            detJ = self.Jacobian(*evaluationPoints)
        return detJ ** 2

    def inverse_Jacobian_matrix(self, *evaluationPoints, J=None):
        """The inverse Jacobian matrix. """
        if J is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        Jacobian = \
            + J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] \
            + J[0][2] * J[1][0] * J[2][1] - J[0][0] * J[1][2] * J[2][1] \
            - J[0][1] * J[1][0] * J[2][2] - J[0][2] * J[1][1] * J[2][0]
        reciprocalJacobian = 1 / Jacobian
        del Jacobian
        iJ00 = reciprocalJacobian * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
        iJ01 = reciprocalJacobian * (J[2][1] * J[0][2] - J[2][2] * J[0][1])
        iJ02 = reciprocalJacobian * (J[0][1] * J[1][2] - J[0][2] * J[1][1])
        iJ10 = reciprocalJacobian * (J[1][2] * J[2][0] - J[1][0] * J[2][2])
        iJ11 = reciprocalJacobian * (J[2][2] * J[0][0] - J[2][0] * J[0][2])
        iJ12 = reciprocalJacobian * (J[0][2] * J[1][0] - J[0][0] * J[1][2])
        iJ20 = reciprocalJacobian * (J[1][0] * J[2][1] - J[1][1] * J[2][0])
        iJ21 = reciprocalJacobian * (J[2][0] * J[0][1] - J[2][1] * J[0][0])
        iJ22 = reciprocalJacobian * (J[0][0] * J[1][1] - J[0][1] * J[1][0])
        return [[iJ00, iJ01, iJ02],
                [iJ10, iJ11, iJ12],
                [iJ20, iJ21, iJ22]]

    def inverse_Jacobian(self, *evaluationPoints, iJ=None):
        """Determinant of the inverse Jacobian matrix. """
        if iJ is None:
            iJ = self.inverse_Jacobian_matrix(*evaluationPoints)
        return + iJ[0][0]*iJ[1][1]*iJ[2][2] \
            + iJ[0][1]*iJ[1][2]*iJ[2][0] + iJ[0][2]*iJ[1][0]*iJ[2][1]\
            - iJ[0][0]*iJ[1][2]*iJ[2][1] - iJ[0][1]*iJ[1][0]*iJ[2][2] - iJ[0][2]*iJ[1][1]*iJ[2][0]

    def metric_matrix(self, *evaluationPoints, J=None):
        """
        Also called metric tensor. Let J be the Jacobian matrix. The ``metricMatrix`` is
        denoted by G, G := J^T.dot(J). And the metric is ``g := (det(J))**2 or g := det(G).``
        Which means for a square Jacobian matrix, the metric turns out to be the square of the
        determinant of the Jacobian matrix.

        The entries of G is normally denoted as g_{i,j}.
        """
        if J is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        G = [[None for _ in range(3)] for __ in range(3)]
        for i in range(3):
            for j in range(i, 3):
                # noinspection PyTypeChecker
                G[i][j] = J[0][i] * J[0][j]
                for L in range(1, 3):
                    G[i][j] += J[L][i] * J[L][j]
                if i != j:
                    G[j][i] = G[i][j]
        return G

    def inverse_metric_matrix(self, *evaluationPoints, iJ=None):
        """
        The ``inverseMetricMatrix`` is the metric matrix of the inverse Jacobian matrix
        or the metric of the inverse mapping. It is usually denoted as G^{-1}.

        The entries of G^{-1} is normally denoted as g^{i,j}.
        """
        if iJ is None:
            iJ = self.inverse_Jacobian_matrix(*evaluationPoints)
        iG = [[None for _ in range(3)] for __ in range(3)]
        for i in range(3):
            for j in range(i, 3):
                # noinspection PyTypeChecker
                iG[i][j] = iJ[i][0] * iJ[j][0]
                for L in range(1, 3):
                    iG[i][j] += iJ[i][L] * iJ[j][L]
                if i != j:
                    iG[j][i] = iG[i][j]
        return iG
