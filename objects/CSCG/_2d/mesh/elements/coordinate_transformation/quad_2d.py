# -*- coding: utf-8 -*-
from components.quadrature import Quadrature
from components.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.elements.coordinate_transformation.helpers.value_cache import ElementsCTValuesCache
import numpy as np
from components.decorators.all import memoize1

class _2dCSCG_ECT_2d_QUAD(FrozenOnly):
    def __init__(self, ect):
        self._elements_ = ect._elements_
        self._freeze_self_()

    @staticmethod
    def ___compute_xietasigma___(quad_degree, quad_type):
        _Quadrature_ = Quadrature(quad_degree, category=quad_type)
        quad_nodes = _Quadrature_.quad[0]
        xi, eta = np.meshgrid(*quad_nodes, indexing='ij')
        return xi, eta

    @memoize1
    def mapping(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta)

    @memoize1
    def X(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta)
    @memoize1
    def Y(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta)

    @memoize1
    def Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta)
    @memoize1
    def J00(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta)
    @memoize1
    def J01(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta)
    @memoize1
    def J10(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta)
    @memoize1
    def J11(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta)

    @memoize1
    def Jacobian(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, intermediateData=None)
    @memoize1
    def metric(self, quad_degree, quad_type):
        """g := det(G) = Jacobian ** 2."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, intermediateData=None)
    @memoize1
    def metric_matrix(self, quad_degree, quad_type):
        """G, g_{i,j}."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, intermediateData=None)

    @memoize1
    def inverse_Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, intermediateData=None)
    @memoize1
    def inverse_Jacobian(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, intermediateData=None)
    @memoize1
    def inverse_metric_matrix(self, quad_degree, quad_type):
        """G^-1, g^{i,j}."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, intermediateData=None)