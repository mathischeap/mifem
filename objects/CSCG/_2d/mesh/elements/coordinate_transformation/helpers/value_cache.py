# -*- coding: utf-8 -*-
"""

"""
from components.freeze.base import FrozenOnly
from root.config.main import CACHE_FACTOR


class ElementsCTValuesCache(FrozenOnly):
    def __init__(self, elements, CTT, xi, eta, intermediateData=None):
        self._elements_ = elements
        self._CTT_ = CTT
        self._xi_eta_sigma_ = (xi, eta)
        if intermediateData is not None:
            assert intermediateData.__class__.__name__ == 'ElementsCTValuesCache', \
                "intermediateData can only be an ElementsCTValues object."
            inter_xi_eta_sigma = intermediateData._xi_eta_sigma_
            for i, com in enumerate(inter_xi_eta_sigma):
                assert self._xi_eta_sigma_[i] is com, \
                    "intermediateData xi, eta sigma must be the same object."
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
        if self._CTT_ in ('mapping', 'X', 'Y'):
            return getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
        elif self._CTT_ in ('Jacobian_matrix', 'J00', 'J01', 'J10', 'J11'):
            if type_wrt_metric in self._cache_:
                return self._cache_[type_wrt_metric]
            else:
                JM = getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
                if isinstance(type_wrt_metric, str) and \
                        type_wrt_metric in self._multi_elements_metric_ and \
                        type_wrt_metric not in self._cache_ and \
                        self._multi_elements_metric_[type_wrt_metric] >= CACHE_FACTOR:
                    self._cache_[type_wrt_metric] = JM
                return JM
        else:
            if type_wrt_metric in self._cache_:
                return self._cache_[type_wrt_metric]
            else:
                if self._intermediateData_ is None:
                    result = getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_)
                else:
                    result = getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_, self._intermediateData_[i])
                if isinstance(type_wrt_metric, str) and \
                        type_wrt_metric in self._multi_elements_metric_ and \
                        type_wrt_metric not in self._cache_ and \
                        self._multi_elements_metric_[type_wrt_metric] >= CACHE_FACTOR:
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
