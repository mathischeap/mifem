

from screws.freeze.base import FrozenOnly
from root.config.main import caChe_factor


class TraceElementsCTValuesCache(FrozenOnly):
    def __init__(self, trace_elements, CTT, evaluation_points_3):
        self._elements_ = trace_elements
        assert isinstance(CTT, str) and CTT != 'mapping', f"CTT={CTT} wrong."
        self._CTT_ = CTT
        self._0ep_ = [evaluation_points_3[1], evaluation_points_3[2]]
        self._e1p_ = [evaluation_points_3[0], evaluation_points_3[2]] # do not use 2, 0
        self._ep2_ = [evaluation_points_3[0], evaluation_points_3[1]]
        self._multi_elements_metric_ = self._elements_._multi_elements_metric_
        self._cache_ = dict()
        self._freeze_self_()

    def __getitem__(self, i):
        element = self._elements_[i]
        type_wrt_metric = element.type_wrt_metric.mark

        if isinstance(type_wrt_metric, int): # it is unique, we use the id (int) as the mark. Otherwise, it must be a string.
            side = element.CHARACTERISTIC_side
            if side in 'NS':
                _xi_eta_sigma_ = self._0ep_
            elif side in 'WE':
                _xi_eta_sigma_ = self._e1p_
            elif side in 'BF':
                _xi_eta_sigma_ = self._ep2_
            else:
                raise Exception()
            return getattr(element.coordinate_transformation, self._CTT_)(
                *_xi_eta_sigma_)

        elif type_wrt_metric in self._cache_:
            return self._cache_[type_wrt_metric]

        else:

            side = element.CHARACTERISTIC_side
            if side in 'NS':
                _xi_eta_sigma_ = self._0ep_
            elif side in 'WE':
                _xi_eta_sigma_ = self._e1p_
            elif side in 'BF':
                _xi_eta_sigma_ = self._ep2_
            else:
                raise Exception()
            result = getattr(element.coordinate_transformation, self._CTT_)(
                *_xi_eta_sigma_)

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


