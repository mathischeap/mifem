# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, the Netherlands

"""


from screws.freeze.base import FrozenOnly



class _2dCSCG_MeshElements_IS(FrozenOnly):
    """"""
    def __init__(self, elements):
        self._elements_ = elements
        self._check_1_ = None
        self._check_2_ = None
        self._freeze_self_()

    @property
    def homogeneous_according_to_types_wrt_metric(self):
        """"""
        if self._check_1_ is None:
            if self._elements_.num <= 1:
                self._check_1_ = True
            else:
                similarity = self._elements_.statistic[
                    'similarity according to types wrt metric']
                self._check_1_ = similarity == 1

        return self._check_1_

    @property
    def all_orthogonal(self):
        """All local elements are orthogonal."""
        if self._check_2_ is None:
            if self._elements_.num < 1:
                self._check_2_ = True
            else:
                num_local_orthogonal_elements = self._elements_.statistic['num_local_orthogonal_elements']
                self._check_2_ = num_local_orthogonal_elements == self._elements_.num

        return self._check_2_