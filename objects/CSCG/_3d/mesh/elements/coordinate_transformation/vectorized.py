# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

import numpy as np
from screws.freeze.base import FrozenOnly


class _3dCSCG_MeshElement_CT_VEC(FrozenOnly):
    """"""
    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()

    def Jacobian(self, xi, et, sg):
        """

        :param xi: ndarray, (xi, et, sg) be same shape.
        :param et: ndarray, (xi, et, sg) be same shape.
        :param sg: ndarray, (xi, et, sg) be same shape.
        :return:
        """
        if len(self._elements_) == 0:
            return None

        if self._elements_.IS.homogeneous_according_to_types_wrt_metric:

            i = self._elements_.indices[0]
            Jacobian = self._elements_[i].coordinate_transformation.Jacobian(xi, et, sg)
            return Jacobian

        else:
            Jacobian = list()

            J = self._elements_.coordinate_transformation.Jacobian(xi, et, sg)

            for i in self._elements_:
                Jacobian.append(J[i])
            return np.array(Jacobian)

    def inverse_Jacobian(self, xi, et, sg):
        """

        :param xi: ndarray, (xi, et, sg) be same shape
        :param et: ndarray, (xi, et, sg) be same shape
        :param sg: ndarray, (xi, et, sg) be same shape
        :return:
        """
        if len(self._elements_) == 0:
            return None

        if self._elements_.IS.homogeneous_according_to_types_wrt_metric:

            i = self._elements_.indices[0]
            inverse_Jacobian = self._elements_[i].coordinate_transformation.inverse_Jacobian(xi, et, sg)
            return inverse_Jacobian

        else:
            inverse_Jacobian = list()

            iJ = self._elements_.coordinate_transformation.inverse_Jacobian(xi, et, sg)

            for i in self._elements_:
                inverse_Jacobian.append(iJ[i])
            return np.array(inverse_Jacobian)