# -*- coding: utf-8 -*-
import sys
if './' not in sys.path:
    sys.path.append('./')

from objects.CSCG._3d.mesh.elements.coordinate_transformation.quad_1d import _3dCSCG_ECT_1d_QUAD
from objects.CSCG._3d.mesh.elements.coordinate_transformation.quad_3d import _3dCSCG_ECT_3d_QUAD
from objects.CSCG._3d.mesh.elements.coordinate_transformation.helpers.value_cache import ElementsCTValuesCache
from objects.CSCG._3d.mesh.elements.coordinate_transformation.vectorized import _3dCSCG_MeshElement_CT_VEC
from components.freeze.main import FrozenOnly


class _3dCSCG_Mesh_Elements_CT(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._vectorized_ = None
        self._ctq_1d_ = _3dCSCG_ECT_1d_QUAD(self)
        self._ctq_3d_ = _3dCSCG_ECT_3d_QUAD(self)
        self._freeze_self_()

    @property
    def vectorized(self):
        if self._vectorized_ is None:
            self._vectorized_ = _3dCSCG_MeshElement_CT_VEC(self._elements_)
        return self._vectorized_

    @property
    def QUAD_1d(self):
        """Evaluating the coordinate transformation from quadrature. Results are in 1d array (do ravel)."""
        return self._ctq_1d_

    @property
    def QUAD_3d(self):
        """Evaluating the coordinate transformation from quadrature. Results are in 3d array."""
        return self._ctq_3d_

    def mapping(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta, sigma)

    def X(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta, sigma)

    def Y(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta, sigma)

    def Z(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'Z', xi, eta, sigma)

    def Jacobian_matrix(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta, sigma)

    def J00(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta, sigma)

    def J01(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta, sigma)

    def J02(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J02', xi, eta, sigma)

    def J10(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta, sigma)

    def J11(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta, sigma)

    def J12(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J12', xi, eta, sigma)

    def J20(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J20', xi, eta, sigma)

    def J21(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J21', xi, eta, sigma)

    def J22(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J22', xi, eta, sigma)

    def Jacobian(self, xi, eta, sigma, J=None):
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, sigma, intermediateData=J)

    def metric(self, xi, eta, sigma, detJ=None):
        """g := det(G) = Jacobian ** 2."""
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, sigma, intermediateData=detJ)

    def metric_matrix(self, xi, eta, sigma, J=None):
        """G, g_{i,j}."""
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, sigma, intermediateData=J)

    def inverse_Jacobian_matrix(self, xi, eta, sigma, J=None):
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, sigma, intermediateData=J)

    def inverse_Jacobian(self, xi, eta, sigma, iJ=None):
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, sigma, intermediateData=iJ)

    def inverse_metric_matrix(self, xi, eta, sigma, iJ=None):
        """G^-1, g^{i,j}."""
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, sigma, intermediateData=iJ)
