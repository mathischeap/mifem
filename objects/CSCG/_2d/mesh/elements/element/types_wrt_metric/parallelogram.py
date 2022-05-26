# -*- coding: utf-8 -*-

import numpy as np
from objects.CSCG._2d.mesh.elements.element.types_wrt_metric.base import ElementTypeWr2MetricBase
from objects.mpRfT._2d.mesh.cell.types_wrt_metric.parallelogram import mpRfT2_ParallelogramCell

from screws.decorators.accepts import accepts


class ParallelogramElement(ElementTypeWr2MetricBase):
    """
    A parallelogramElement element is an element that:
        1) four edges are straight lines and form a parallelogram;
        2) internal transformation is linear,
        3) it is not an OrthogonalElement element.

    So, we know, it can be a rectangle or square once its left edge is not parallel with x-axis. If
    its left edge is parallel with x-axis, such a rectangle or square is an orthogonal element.
    """
    @accepts('self', (float, int), (float, int), (float, int), (float, int))
    def __init__(self, angleL, L, L_angle_U, U):
        self._data_ = angleL, L, L_angle_U, U
        assert 0 < L_angle_U < np.pi, f"L.U angle = {L_angle_U} wrong!"
        assert 0 <= angleL < 2 * np.pi, f" angle L={angleL} must be < 2pi."
        self._mark_ = 'Parallelogram.' + 'aL{}_L{}_{}_U{}'.format(
            '%.5f' % angleL, '%.5f' % L, '%.5f' % L_angle_U, '%.5f' % U)
        self._freeze_self_()


    def ___CLASSIFY_mpRfT2_CELL_of_origin_and_delta___(self, origin_and_delta):
        """"""
        angleL, L, L_angle_U, U = self._data_
        delta = origin_and_delta[1]
        L *= delta / 2
        U *= delta / 2
        return mpRfT2_ParallelogramCell(angleL, L, L_angle_U, U)