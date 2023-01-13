# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from abc import ABC

from objects.CSCG._2d.fields.base import _2dCSCG_Continuous_FORM_BASE


class _2dCSCG_TensorField(_2dCSCG_Continuous_FORM_BASE, ABC, ndim=2):
    """The continuous vector field."""

    def __repr__(self):
        """"""
        return f"2dCSCG_tensor_field=<{self.standard_properties.name}>@{id(self)}"

    @property
    def shape(self):
        return 2, 2
