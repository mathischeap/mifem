# -*- coding: utf-8 -*-
"""
The type classes for trace elements w.r.t. metric.

"""

from screws.freeze.main import FrozenOnly






class TraceElementTypeWr2MetricBase(FrozenOnly):
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