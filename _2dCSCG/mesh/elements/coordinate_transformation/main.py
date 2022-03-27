

from screws.freeze.main import FrozenOnly
from _2dCSCG.mesh.elements.coordinate_transformation.helpers.value_cache import ElementsCTValuesCache
from _2dCSCG.mesh.elements.coordinate_transformation.quad_1d import _2dCSCG_ECT_1d_QUAD
from _2dCSCG.mesh.elements.coordinate_transformation.quad_2d import _2dCSCG_ECT_2d_QUAD
from _2dCSCG.mesh.elements.coordinate_transformation.vectorized import _2dCSCG_MeshElements_CT_Vectorized


class _2dCSCG_Mesh_Elements_CT(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._vectorized_ = None
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._ctq_1d_ = _2dCSCG_ECT_1d_QUAD(self)
        self._ctq_2d_ = _2dCSCG_ECT_2d_QUAD(self)

    @property
    def vectorized(self):
        if self._vectorized_ is None:
            self._vectorized_ = _2dCSCG_MeshElements_CT_Vectorized(self._elements_)
        return self._vectorized_

    @property
    def QUAD_1d(self):
        """Evaluating the coordinate transformation from quadrature. Results are in 1d array (do ravel)."""
        return self._ctq_1d_

    @property
    def QUAD_2d(self):
        """Evaluating the coordinate transformation from quadrature. Results are in 2d array."""
        return self._ctq_2d_


    def mapping(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta)
    def X(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta)
    def Y(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta)

    def Jacobian_matrix(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta)
    def J00(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta)
    def J01(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta)
    def J10(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta)
    def J11(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta)

    def Jacobian(self, xi, eta, J=None):
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, intermediateData=J)
    def metric(self, xi, eta, detJ=None):
        """g := det(G) = Jacobian ** 2."""
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, intermediateData=detJ)
    def metric_matrix(self, xi, eta, J=None):
        """G, g_{i,j}."""
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, intermediateData=J)

    def inverse_Jacobian_matrix(self, xi, eta, J=None):
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, intermediateData=J)
    def inverse_Jacobian(self, xi, eta, iJ=None):
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, intermediateData=iJ)
    def inverse_metric_matrix(self, xi, eta, iJ=None):
        """G^-1, g^{i,j}."""
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, intermediateData=iJ)
