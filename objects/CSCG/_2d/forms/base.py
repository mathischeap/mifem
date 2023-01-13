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
        self.___kwargs___ = None
        assert self.ndim == 2, "CHECK ndim"

    @property
    def dqp(self):
        """Return the default quadrature degree for variant integrations."""
        if self.space.__class__.__name__ == '_2dCSCG_PolynomialSpace':
            return self.p
        else:
            raise NotImplementedError()

    def __add__(self, other):
        """"""
        assert other.__class__.__name__ == self.__class__.__name__, f"forms do not match."
        assert other.mesh == self.mesh, f"meshes do not match."
        assert other.space == self.space, f"spaces do not match."

        assert self.cochain.local is not None,  f"a of (a+b) has no cochain.local, cannot perform add operator."
        assert other.cochain.local is not None, f"b of (a+b) has no cochain.local, cannot perform add operator."

        kwargs_A = self.___kwargs___
        kwargs_B = other.___kwargs___

        if 'name' in kwargs_A:
            name_A = kwargs_A['name']
        else:
            name_A = 'form_A'
        if 'name' in kwargs_B:
            name_B = kwargs_B['name']
        else:
            name_B = 'form_B'

        name = name_A + "+" + name_B

        kwargs = dict()
        for key in kwargs_A:
            if key == 'name':
                kwargs['name'] = name
            else:
                kwargs[key] = kwargs_A[key]

        if 'name' not in kwargs:
            kwargs['name'] = name

        # noinspection PyArgumentList
        result_form = self.__class__(self.mesh, self.space, **kwargs)

        # MUST do: add the cochain local ------------------
        COCHAIN_LOCAL = dict()
        for e in self.mesh.elements:
            clA = self.cochain.local[e]
            clB = other.cochain.local[e]
            COCHAIN_LOCAL[e] = clA + clB

        # ========================================================
        result_form.cochain.local = COCHAIN_LOCAL
        return result_form

    def __rmul__(self, other):
        """other * self

        Parameters
        ----------
        other

        Returns
        -------

        """

        if isinstance(other, (int, float)):

            kwargs = dict()
            for key in self.___kwargs___:
                kwargs[key] = self.___kwargs___[key]

            if 'name' in kwargs:
                kwargs['name'] = f"{other}*{kwargs['name']}"
            else:
                kwargs['name'] = f"{other}*form"

            result_form = self.__class__(self.mesh, self.space, **kwargs)

            COCHAIN_LOCAL = dict()
            for e in self.mesh.elements:
                clA = self.cochain.local[e]
                COCHAIN_LOCAL[e] = other * clA

            # ========================================================
            result_form.cochain.local = COCHAIN_LOCAL
            return result_form

        else:
            raise NotImplementedError()
