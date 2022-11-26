# -*- coding: utf-8 -*-


from components.freeze.main import FrozenOnly


class ElementTypeWr2MetricBase(FrozenOnly):
    """A base for all element types w.r.t. metric.

    For each type of regions, we can classify its elements into different types. These types are all coded here.
    """
    @property
    def mark(self):
        """A mark is key that identifies the element metric.

        If the marks of two elements are the same, then they have the same metric, otherwise,
        their metric are different. A mark normally is a string. But for chaotic element,
        it is an int: the id of the object.

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

    def ___CLASSIFY_3nCSCG_RF2_CELL_of_origin_and_delta___(self, origin_and_delta):
        raise NotImplementedError(
            f"Please implement ___CLASSIFY_3nCSCG_RF2_CELL_of_origin_and_delta___ for "
            f"TypeWr2Metric named: {self.__class__.__name__}"
        )