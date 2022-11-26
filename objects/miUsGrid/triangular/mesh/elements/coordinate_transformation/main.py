# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/6/2022 6:03 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly

from objects.miUsGrid.triangular.mesh.elements.coordinate_transformation.helpers.cache import \
    miUsTriangle_Elements_CT_Cache


class miUsTriangle_Elements_CoordinateTransformation(FrozenOnly):
    """"""

    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()


    def mapping(self, r, s):
        return miUsTriangle_Elements_CT_Cache(
            self._elements_,
            'mapping',
            r, s,
        )

    def Jacobian_matrix(self, r, s):
        return miUsTriangle_Elements_CT_Cache(
            self._elements_,
            'Jacobian_matrix',
            r, s,
        )

    def Jacobian(self, r, s, itmD=None):
        """Determinant of the Jacobian matrix.

        intermediate_data (itmD) is: Jacobian_matrix
        """
        return miUsTriangle_Elements_CT_Cache(
            self._elements_,
            'Jacobian',
            r, s,
            intermediate_data=itmD
        )

    def metric(self, r, s, itmD=None):
        """
        The metric ``g:= det(G):=(det(J))**2``. Since our Jacobian and inverse of Jacobian are both square,
        we know that the metric ``g`` is equal to square of ``det(J)``. ``g = (det(J))**2`` is due to the
        fact that the Jacobian matrix is square. The definition of ``g`` usually is given
        as ``g:= det(G)`` where ``G`` is the metric matrix, or metric tensor.

        intermediate_data (itmD) is: Jacobian
        """
        return miUsTriangle_Elements_CT_Cache(
            self._elements_,
            'metric',
            r, s,
            intermediate_data=itmD
        )

    def inverse_Jacobian_matrix(self, r, s, itmD=None):
        """The inverse Jacobian matrix.

        intermediate_data (itmD) is: Jacobian_matrix
        """
        return miUsTriangle_Elements_CT_Cache(
            self._elements_,
            'inverse_Jacobian_matrix',
            r, s,
            intermediate_data=itmD
        )

    def inverse_Jacobian(self, r, s, itmD=None):
        """Determinant of the inverse Jacobian matrix.

        intermediate_data (itmD) is: inverse_Jacobian_matrix
        """
        return miUsTriangle_Elements_CT_Cache(
            self._elements_,
            'inverse_Jacobian',
            r, s,
            intermediate_data=itmD
        )

    def metric_matrix(self, r, s, itmD=None):
        """
        Also called metric tensor. Let J be the Jacobian matrix. The ``metricMatrix`` is
        denoted by G, G := J^T.dot(J). And the metric is ``g := (det(J))**2 or g := det(G).``
        Which means for a square Jacobian matrix, the metric turns out to be the square of the
        determinant of the Jacobian matrix.

        The entries of G are normally denoted as g_{i,j}.

        intermediate_data (itmD) is: Jacobian_matrix
        """
        return miUsTriangle_Elements_CT_Cache(
            self._elements_,
            'metric_matrix',
            r, s,
            intermediate_data=itmD
        )

    def inverse_metric_matrix(self, r, s, itmD=None):
        """
        The ``inverseMetricMatrix`` is the metric matrix of the inverse Jacobian matrix
        or the metric of the inverse mapping. It is usually denoted as G^{-1}.

        The entries of G^{-1} are normally denoted as g^{i,j}.

        intermediate_data (itmD) is: inverse_Jacobian_matrix
        """
        return miUsTriangle_Elements_CT_Cache(
            self._elements_,
            'inverse_metric_matrix',
            r, s,
            intermediate_data=itmD
        )



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/elements/coordinate_transformation/main.py
    pass
