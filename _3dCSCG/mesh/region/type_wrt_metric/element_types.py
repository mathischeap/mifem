# -*- coding: utf-8 -*-


from SCREWS.frozen import FrozenOnly
from SCREWS.decorators import accepts


class ElementTypeWr2Metric(FrozenOnly):
    """
    A base for all element types w.r.t. metric. For each type of region, we can classify its elements
    into different types. These types are all coded here.
    """
    @property
    def mark(self):
        """
        A mark is key that identifies the element metric. If the marks of two elements are the same, then
        they have the same metric, otherwise, their metric are different. A mark normally is a string. But for
        chaotic element, it is an int: the id of the object.

        :return:
        """
        # noinspection PyUnresolvedReferences
        return self._mark_

    def __eq__(self, other):
        """We ask that the marks to be the same."""
        # The later judge is to make sure we are not comparing it to something else having a mark property
        return self.mark == other.mark and other.___IS_3dCSCG_ElementTypeWr2Metric___

    @property
    def ___IS_3dCSCG_ElementTypeWr2Metric___(self):
        return True


# -------- particular element types below --------------------------------------------------------------------


class ChaoticElement(ElementTypeWr2Metric):
    """
    Chaotic element is the element that its metric is unique. So no any other elements can be the same with
    it. To make sure that happens, we use the unique id as its mark.
    """
    def __init__(self):
        self._mark_ = id(self)
        self._freeze_self_()




class OrthogonalElement(ElementTypeWr2Metric):
    """
    For orthogonal element, we only need to know its lengths on three directions, then we can fix its metric.
    """
    @accepts('self', (tuple, list, 'ndarray', 'ndim=1', 'shape=(3)'))
    def __init__(self, LxLyLz):
        Lx, Ly, Lz = LxLyLz
        self._mark_ = 'Orth.' + 'x{}y{}z{}'.format('%.6f' % Lx, '%.6f' % Ly, '%.6f' % Lz)
        self._freeze_self_()
