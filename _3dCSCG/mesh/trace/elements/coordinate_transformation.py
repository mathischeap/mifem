
import sys
if './' not in sys.path: sys.path.append('./')

from SCREWS.frozen import FrozenOnly
from root.config import caChe_factor


class _3dCSCG_Trace_Elements_CoordinateTransformation(FrozenOnly):
    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def mapping(self, *args, **kwargs):
        """"""
        raise Exception("Can not access mapping from trace elements, must do it from particular trace element.")

    def Jacobian_matrix(self, *evaluation_points_3):
        """
        The local Jacobian matrix.

        :param evaluation_points_3: A tuple or list of shape (ndim, ...).
        """
        return TraceElementsCTValuesCache(self._elements_, 'Jacobian_matrix', evaluation_points_3)

    def inverse_Jacobian_matrix(self, *evaluation_points_3):
        """
        The local inverse_Jacobian matrix.

        :param evaluation_points_3 : A tuple or list of shape (ndim, ...).
        """
        return TraceElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', evaluation_points_3)

    def metric_matrix(self, *evaluation_points_3):
        """Compute the metric matrix G."""
        return TraceElementsCTValuesCache(self._elements_, 'metric_matrix', evaluation_points_3)

    def metric(self, *evaluation_points_3):
        """return metric g."""
        return TraceElementsCTValuesCache(self._elements_, 'metric', evaluation_points_3)

    def unit_normal_vector(self, *evaluation_points_3):
        """return metric g."""
        return TraceElementsCTValuesCache(self._elements_, 'unit_normal_vector', evaluation_points_3)

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






if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace\elements\coordinate_transformation.py
    from _3dCSCG.main import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.SELFCHECK.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality

    mesh.trace.elements.DO.illustrate_trace_element(1)