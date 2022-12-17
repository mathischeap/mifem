# -*- coding: utf-8 -*-
from objects.CSCG._2d.mesh.trace.elements.element.types_wrt_metric.base import TraceElementTypeWr2MetricBase


class OrthogonalTraceElement(TraceElementTypeWr2MetricBase):
    """
    An orthogonal trace element must be a rectangle or square, and it is
    perpendicular to one of the three axes. And has no rotation. That means each edge must be
    perpendicular to the corresponding axes.

    For an orthogonal trace element, we only need to know the axis it is
    perpendicular to and its size, then we can fix its metric.
    """
    def __init__(self, along, length):
        """
        :param along: 'x', 'y'
        :param length: The length of this 1d orthogonal trace element.
        """
        assert len(along) == 1 and along in 'xy', f"along={along} wrong."
        assert isinstance(length, (int, float)) and length > 0, f"length={length} wrong."
        self._mark_ = f'Orth.{along}' + '{}'.format('%.8f' % length)
        self._freeze_self_()
