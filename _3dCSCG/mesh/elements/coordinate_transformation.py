



import sys
if './' not in sys.path: sys.path.append('./')


from root.config.main import *
from screws.freeze.main import FrozenOnly
from screws.quadrature import Quadrature



class _3dCSCG_Mesh_Elements_CT(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._ctq_1d_ = _3dCSCG_ECT_1d_QUAD(self)
        self._ctq_3d_ = _3dCSCG_ECT_3d_QUAD(self)

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


class ElementsCTValuesCache(FrozenOnly):
    def __init__(self, elements, CTT, xi, eta, sigma, intermediateData=None):
        self._elements_ = elements
        self._CTT_ = CTT
        self._xi_eta_sigma_ = (xi, eta, sigma)
        if intermediateData is not None:
            assert intermediateData.__class__.__name__ == 'ElementsCTValuesCache', \
                "intermediateData can only be an ElementsCTValues object."
            inter_xi_eta_sigma = intermediateData._xi_eta_sigma_
            for i, com in enumerate(inter_xi_eta_sigma):
                assert self._xi_eta_sigma_[i] is com, \
                    "intermediateData xi, eta, sigma must be the same object."
            if CTT == 'Jacobian':
                assert intermediateData._CTT_ == 'Jacobian_matrix'
            elif CTT == 'metric':
                assert intermediateData._CTT_ == 'Jacobian'
            elif CTT == 'metric_matrix':
                assert intermediateData._CTT_ == 'Jacobian_matrix'
            elif CTT == 'inverse_Jacobian_matrix':
                assert intermediateData._CTT_ == 'Jacobian_matrix'
            elif CTT == 'inverse_Jacobian':
                assert intermediateData._CTT_ == 'inverse_Jacobian_matrix'
            elif CTT == 'inverse_metric_matrix':
                assert intermediateData._CTT_ == 'inverse_Jacobian_matrix'
            else:
                raise Exception(f"{CTT} ask for no intermediateData or provided wrong intermediateData.")
        else:
            pass
        self._intermediateData_ = intermediateData
        self._multi_elements_metric_ = self._elements_._multi_elements_metric_
        self._cache_ = dict()
        self._freeze_self_()

    def __getitem__(self, i):
        element = self._elements_[i]
        type_wrt_metric = element.type_wrt_metric.mark
        if self._CTT_ in ('mapping', 'X', 'Y', 'Z'):
            return getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
        elif self._CTT_ in ('Jacobian_matrix', 'J00', 'J01', 'J02', 'J10', 'J11', 'J12', 'J20', 'J21', 'J22'):

            if isinstance(type_wrt_metric, int):
                return getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
            elif type_wrt_metric in self._cache_:
                return self._cache_[type_wrt_metric]
            else:
                JM = getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
                if type_wrt_metric in self._multi_elements_metric_ and \
                    type_wrt_metric not in self._cache_ and \
                    self._multi_elements_metric_[type_wrt_metric] >= caChe_factor:
                    self._cache_[type_wrt_metric] = JM
                return JM
        else:
            if isinstance(type_wrt_metric, int):
                if self._intermediateData_ is None:
                    return getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_)
                else:
                    return getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_, self._intermediateData_[i])
            elif type_wrt_metric in self._cache_:
                return self._cache_[type_wrt_metric]
            else:
                if self._intermediateData_ is None:
                    result = getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_)
                else:
                    result = getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_, self._intermediateData_[i])
                if type_wrt_metric in self._multi_elements_metric_ and \
                    type_wrt_metric not in self._cache_ and \
                    self._multi_elements_metric_[type_wrt_metric] >= caChe_factor:
                    # here we have very strict cache rule.
                    self._cache_[type_wrt_metric] = result
                return result

    def __len__(self):
        return len(self._elements_)

    def __contains__(self, item):
        return item in self._elements_

    def __iter__(self):
        for i in self._elements_:
            yield i


class _3dCSCG_ECT_1d_QUAD(FrozenOnly):
    def __init__(self, ect):
        self._elements_ = ect._elements_
        self._freeze_self_()

    @staticmethod
    def ___compute_xietasigma___(quad_degree, quad_type):
        _Quadrature_ = Quadrature(quad_degree, category=quad_type)
        quad_nodes = _Quadrature_.quad[0]
        xi, eta, sigma = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        sigma = sigma.ravel('F')
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










