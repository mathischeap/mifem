
import numpy as np
from screws.quadrature import Quadrature
from screws.freeze.base import FrozenOnly
from _3dCSCG.mesh.elements.coordinate_transformation.helpers.value_cache import ElementsCTValuesCache

class _3dCSCG_ECT_3d_QUAD(FrozenOnly):
    def __init__(self, ect):
        self._elements_ = ect._elements_
        self._freeze_self_()

    @staticmethod
    def ___compute_xietasigma___(quad_degree, quad_type):
        _Quadrature_ = Quadrature(quad_degree, category=quad_type)
        quad_nodes = _Quadrature_.quad[0]
        xi, eta, sigma = np.meshgrid(*quad_nodes, indexing='ij')
        return xi, eta, sigma

    def mapping(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta, sigma)

    def X(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta, sigma)
    def Y(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta, sigma)
    def Z(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Z', xi, eta, sigma)

    def Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta, sigma)
    def J00(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta, sigma)
    def J01(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta, sigma)
    def J02(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J02', xi, eta, sigma)
    def J10(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta, sigma)
    def J11(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta, sigma)
    def J12(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J12', xi, eta, sigma)
    def J20(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J20', xi, eta, sigma)
    def J21(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J21', xi, eta, sigma)
    def J22(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J22', xi, eta, sigma)

    def Jacobian(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, sigma, intermediateData=None)
    def metric(self, quad_degree, quad_type):
        """g := det(G) = Jacobian ** 2."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, sigma, intermediateData=None)
    def metric_matrix(self, quad_degree, quad_type):
        """G, g_{i,j}."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, sigma, intermediateData=None)

    def inverse_Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, sigma, intermediateData=None)
    def inverse_Jacobian(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, sigma, intermediateData=None)
    def inverse_metric_matrix(self, quad_degree, quad_type):
        """G^-1, g^{i,j}."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, sigma, intermediateData=None)


