



import numpy as np
from _2dCSCG.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase


from screws.decorators import accepts


class ParallelogramElement(ElementTypeWr2MetricBase):
    """
    A parallelogramElement element is an element that: 1) four edges are straight line; 2) internal
    transformation is linear, 3) it is not an OrthogonalElement element.

    So, we know, it can be a rectangle or square once its left edge is not parallel with x-axis. If
    its left edge is parallel with x-axis, such a rectangle or square is an orthogonal element.
    """
    @accepts('self', (float, int), (float, int), (float, int), (float, int))
    def __init__(self, angleL, L, L_angle_U, U):
        assert 0 < L_angle_U < np.pi, f"L.U angle = {L_angle_U} wrong!"
        assert 0 <= angleL < 2 * np.pi, f" angle L={angleL} must be < 2pi."
        self._mark_ = 'Parallelogram.' + 'aL{}_L{}_{}_U{}'.format(
            '%.5f' % angleL, '%.5f' % L, '%.5f' % L_angle_U, '%.5f' % U)
        self._freeze_self_()