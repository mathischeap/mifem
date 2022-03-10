# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.freeze.main import FrozenOnly



class TypeWr2Metric(FrozenOnly):
    def __init__(self, region):
        self._region_ = region
        self._mark_ = None

    def __eq__(self, other):
        return self.mark == other.mark and other.___IS_2dCSCG_TypeWr2Metric___
        # make sure we can comparing regions type metric objects.

    @property
    def ___IS_2dCSCG_TypeWr2Metric___(self):
        return True

    @property
    def mark(self):
        raise NotImplementedError(
            f"Please implement property mark for TypeWr2Metric named: {self.__class__.__name__}"
        )

    def ___CLASSIFY_ELEMENT_of_spacing___(self, spacing: tuple):
        raise NotImplementedError(
            f"Please implement ___CLASSIFY_ELEMENT_of_spacing___ for "
            f"TypeWr2Metric named: {self.__class__.__name__}"
        )




# -----------------------------------------------------------------------------------------------------------



