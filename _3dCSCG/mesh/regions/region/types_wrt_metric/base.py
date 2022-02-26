# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.frozen import FrozenOnly



class TypeWr2MetricBase(FrozenOnly):
    """A parent of all regions types w.r.t. metric."""
    def __init__(self, region):
        self._region_ = region
        self._mark_ = None

    def __eq__(self, other):
        return self.mark == other.mark and other.___IS_3dCSCG_TypeWr2Metric___
        # make sure we can compare regions type metric objects.

    @property
    def ___IS_3dCSCG_TypeWr2Metric___(self):
        return True

    @property
    def mark(self):
        raise NotImplementedError(
            f"Please implement property mark for TypeWr2Metric named: {self.__class__.__name__}"
        )

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple):
        """

        :param spacing: The spacing that representing the element.
        :return:
        """
        raise NotImplementedError(
            f"Please implement ___CLASSIFY_ELEMENT_of_spacing___ for "
            f"TypeWr2Metric named: {self.__class__.__name__}"
        )

    def ___CLASSIFY_TRACE_ELEMENT_of_spacing___(self, trace_spacing: tuple):
        """

        :param trace_spacing: the trace_spacing representing a trace element.
        :return:
        """
        raise NotImplementedError(
            f"Please implement ___CLASSIFY_TRACE_ELEMENT_of_spacing___ for "
            f"TypeWr2Metric named: {self.__class__.__name__}"
        )

