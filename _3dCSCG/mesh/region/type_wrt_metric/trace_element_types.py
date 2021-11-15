# -*- coding: utf-8 -*-
"""
The type classes for trace elements w.r.t. metric.

"""

from SCREWS.frozen import FrozenOnly






class TraceElementTypeWr2Metric(FrozenOnly):
    """
    A base for all trace element types w.r.t. metric. For each type of mesh element, we can classify its sides (trace elements)
    into different types. These types are all coded here.
    """
    @property
    def mark(self):
        """
        A mark is key that identifies the trace element metric. If the marks of two trace elements are the same, then
        they have the same metric, otherwise, their metric are different. A mark normally is a string. But for
        chaotic trace element, it is an int: the id of the object.

        :return:
        """
        # noinspection PyUnresolvedReferences
        return self._mark_

    def __eq__(self, other):
        """We ask that the marks to be the same."""
        # The later judge is to make sure we are not comparing it to something else having a mark property
        return self.mark == other.mark and other.___IS_3dCSCG_TraceElementTypeWr2Metric___

    @property
    def ___IS_3dCSCG_TraceElementTypeWr2Metric___(self):
        return True




# -------- particular trace element types below --------------------------------------------------------------------


class ChaoticTraceElement(TraceElementTypeWr2Metric):
    """
    Chaotic trace element is the trace element that its metric is unique. So no any other trace elements can be the same with
    it. To make sure that happens, we use the unique id as its mark.
    """
    def __init__(self):
        self._mark_ = id(self)
        self._freeze_self_()




class OrthogonalTraceElement(TraceElementTypeWr2Metric):
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

