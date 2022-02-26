# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from _2dCSCG.fields.base import _2dCSCG_Continuous_FORM_BASE


class _2dCSCG_TensorField(_2dCSCG_Continuous_FORM_BASE, ndim=2):
    """The continuous vector field."""

    @property
    def shape(self):
        return (2, 2)

