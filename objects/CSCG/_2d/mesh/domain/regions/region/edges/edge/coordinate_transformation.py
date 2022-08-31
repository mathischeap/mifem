# -*- coding: utf-8 -*-

import numpy as np
from screws.freeze.main import FrozenOnly


class EdgeCoordinateTransformation(FrozenOnly):
    """[0,1] -> edge """

    def __init__(self, edge):
        """ """
        assert edge.ndim == 2, " <Region> <DII> <EdgeCoordinateTransformation> "
        assert edge.__class__.__name__ == 'Edge', " <Region> <DII> <EdgeCoordinateTransformation> "
        self._edge_ = edge
        self._freeze_self_()

    @property
    def ndim(self):
        return 2

    @staticmethod
    def ___check_t___(t):
        """ """
        assert np.ndim(t) == 1, " <EdgeCoordinateTransformation> : 0 <= t <= 1."
        assert np.min(t) >= 0 and np.max(t) <= 1, " <EdgeCoordinateTransformation> : 0 <= t <= 1."
        if np.size(t) > 1:
            assert all(np.diff(t) > 0), " <EdgeCoordinateTransformation> : 0 <= t <= 1."
        return t


    def mapping(self, t):
        """
        Parameters
        ----------
        t :
            0 <= t <= 1.

        """
        t = self.___check_t___(t)
        O = np.zeros(np.shape(t))
        I = np.ones(np.shape(t))

        if self._edge_._name_ == 'U':  # U
            return self._edge_._region_.interpolation.mapping(O, t)
        elif self._edge_._name_ == 'D':  # S
            return self._edge_._region_.interpolation.mapping(I, t)
        elif self._edge_._name_ == 'L':  # L
            return self._edge_._region_.interpolation.mapping(t, O)
        elif self._edge_._name_ == 'R':  # R
            return self._edge_._region_.interpolation.mapping(t, I)
        else:
            raise Exception()

    def tangent_vector(self, t):
        """
        It is the Jacobian matrix as well.

        Parameters
        ----------
        t :
            0 <= t <= 1.

        """
        t = self.___check_t___(t)
        O = np.zeros(np.shape(t))
        I = np.ones(np.shape(t))
        if self._edge_._name_ == 'U':  # U
            FJM = self._edge_._region_.interpolation.Jacobian_matrix(O, t)
            return FJM[0][1], FJM[1][1]
        elif self._edge_._name_ == 'D':  # S
            FJM = self._edge_._region_.interpolation.Jacobian_matrix(I, t)
            return FJM[0][1], FJM[1][1]
        elif self._edge_._name_ == 'L':  # L
            FJM = self._edge_._region_.interpolation.Jacobian_matrix(t, O)
            return FJM[0][0], FJM[1][0]
        elif self._edge_._name_ == 'R':  # R
            FJM = self._edge_._region_.interpolation.Jacobian_matrix(t, I)
            return FJM[0][0], FJM[1][0]
        else:
            raise Exception()

    def unit_tangent_vector(self, t):
        """ """
        Jm = self.tangent_vector(t)
        norm = np.sqrt(Jm[0] ** 2 + Jm[1] ** 2)
        return Jm[0] / norm, Jm[1] / norm

    def normal_vector(self, t):
        """
        Note that this is not always an out-ward normal vector. On Right and Upper
        edges, it is. But, on Left and Down, it in in-ward.

        """
        tv = self.tangent_vector(t)
        rotated_tv = [-tv[1], tv[0]]  # we always rotate anti_clock_wise 90 degree.
        return rotated_tv

    def unit_normal_vector(self, t):
        """
        Note that this is not always an unit out-ward normal vector. On Right and Upper
        edges, it is. But, on Left and Down, it in in-ward.

        """
        nv = self.normal_vector(t)
        norm = np.sqrt(nv[0] ** 2 + nv[1] ** 2)
        return nv[0] / norm, nv[1] / norm