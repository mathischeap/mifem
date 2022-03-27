
import sys
if './' not in sys.path: sys.path.append('../')

from screws.freeze.main import FrozenOnly
from _3dCSCG.mesh.trace.elements.coordinate_transformation.helpers.value_cache import TraceElementsCTValuesCache


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




if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace\elements\coordinate_transformation.py
    from _3dCSCG.master import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.SELFCHECK.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality

    mesh.trace.elements.do.illustrate_trace_element(1)