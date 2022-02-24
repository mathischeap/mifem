

from _3dCSCG.mesh.trace.elements.element.types_wrt_metric.base import TraceElementTypeWr2MetricBase




class OrthogonalTraceElement(TraceElementTypeWr2MetricBase):
    """
    An orthogonal trace element must be a rectangle or square and it is
    perpendicular to one of the three axes. And has no rotation. That means each edge must be
    perpendicular to the corresponding axes.

    For an orthogonal trace element, we only need to know the axis it is
    perpendicular to and its size, then we can fix its metric.
    """
    def __init__(self, perp_to, d1, d2):
        """
        :param perp_to: 'x', 'y' or 'z'
        :param d1: if perp_to='x', d1=dy, d2=dz, else if perp_to='y', d1
            =dz, d2=dx, else (perp_to='z'), d1=dx, d2=dy.
        :param d2:
        """
        assert len(perp_to) == 1 and perp_to in 'xyz', f"perp_to={perp_to} wrong."
        assert isinstance(d1, (int, float)) and d1 > 0, f"d1={d1} wrong."
        assert isinstance(d2, (int, float)) and d2 > 0, f"d2={d2} wrong."
        self._mark_ = f'Orth.{perp_to}' + 'd1{}d2{}'.format('%.6f' % d1, '%.6f' % d2)
        self._freeze_self_()
