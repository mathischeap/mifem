# -*- coding: utf-8 -*-
"""
CONTINUOUS FORMS are special, they do not inherit the main class of form!

"""
from abc import ABC
from objects.CSCG.base.fields.base import CSCG_Continuous_FORM_BASE


class _3dCSCG_Continuous_FORM_BASE(CSCG_Continuous_FORM_BASE, ABC):
    """"""
    def __init_subclass__(cls, ndim=3):
        super().__init_subclass__(ndim=ndim)

    def __init__(self, mesh, ftype, valid_time):
        """

        Parameters
        ----------
        mesh
        ftype
        valid_time
        """
        assert mesh.__class__.__name__ == '_3dCSCG_Mesh', "Need a 3dCSCG mesh."
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_field')
        super().__init__(mesh, ftype, valid_time)
