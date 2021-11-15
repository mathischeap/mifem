# -*- coding: utf-8 -*-

import numpy as np
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
        return self.mark == other.mark and other.___IS_2dCSCG_ElementTypeWr2Metric___

    @property
    def ___IS_2dCSCG_ElementTypeWr2Metric___(self):
        return True



# -------- particular element types below --------------------------------------------------------------------


class ChaoticElement(ElementTypeWr2Metric):
    """
    Chaotic element is the element that its metric is unique. So no any other elemenets can be the same with
    it. To make sure that happens, we use the unique id as its mark.
    """
    def __init__(self):
        self._mark_ = id(self)
        self._freeze_self_()




class OrthogonalElement(ElementTypeWr2Metric):
    """
    An orthogonal element is an element: 1) it is a rectangle (including square), 2) its internal
    transformation is linear, 3) its left edge is parallel with x-axis.
    """
    @accepts('self', (tuple, list, 'ndarray', 'ndim=1', 'shape=(2)'))
    def __init__(self, LxLy):
        Lx, Ly = LxLy
        self._mark_ = 'Orth.' + 'x{}y{}'.format('%.3f' % Lx, '%.3f' % Ly)
        self._freeze_self_()





class ParallelogramElement(ElementTypeWr2Metric):
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