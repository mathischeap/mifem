# -*- coding: utf-8 -*-
"""

A BASE for all forms except continuous forms.

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from objects.CSCG.base.forms.base.main import CSCG_FORM_BASE


# noinspection PyAbstractClass
class _2dCSCG_FORM_BASE(CSCG_FORM_BASE):
    """This a parent for all forms. It initializes some fundamental properties,
    like, ``mesh``, ``space``, ``ndim``, ``p`` (basis function degree) and
    ``defaultQuadDegree (dqp)``.
    """
    def __init_subclass__(cls, ndim=2):
        super().__init_subclass__(ndim=ndim)
        cls.___ndim___ = ndim

    def __init__(self, mesh, space, name):
        assert mesh.__class__.__name__ == '_2dCSCG_Mesh', "Need a 2dCSCG mesh."
        assert '2dCSCG|structured|space' in space.standard_properties.stamp, "Need a 2dCSCG space."
        assert mesh.ndim == space.ndim == 2
        super().__init__(mesh, space, name)
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_form')
        assert self.ndim == 2, "CHECK ndim"


    @property
    def dqp(self):
        """Return the default quadrature degree for variant integrations."""
        if self.space.__class__.__name__ == '_2dCSCG_PolynomialSpace':
            return self.p
        else:
            raise NotImplementedError()