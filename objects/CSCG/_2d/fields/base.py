# -*- coding: utf-8 -*-
"""
CONTINUOUS FORMS are special, they do not inherit the main class of form!

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""


from abc import ABC
from objects.CSCG.base.fields.base import CSCG_Continuous_FORM_BASE




class _2dCSCG_Continuous_FORM_BASE(CSCG_Continuous_FORM_BASE, ABC):
    """"""
    def __init_subclass__(cls, ndim=2):
        super().__init_subclass__(ndim=ndim)

    def __init__(self, mesh, ftype, valid_time):
        """"""
        assert mesh.__class__.__name__ == '_2dCSCG_Mesh', "Need a 2dCSCG mesh."
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_field')
        super().__init__(mesh, ftype, valid_time)