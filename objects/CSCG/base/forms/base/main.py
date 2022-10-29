# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from screws.freeze.main import FrozenClass

from objects.CSCG.base.forms.base.BC.main import CSCG_Form_BC



class CSCG_FORM_BASE(FrozenClass):
    """"""
    def __init_subclass__(cls, ndim=None):
        super().__init_subclass__()
        cls.___ndim___ = ndim

    def __init__(self, mesh, space):
        self._mesh_ = mesh
        self._space_ = space
        self._k_ = None
        self.___define_parameters___ = None
        self.standard_properties.___PRIVATE_add_tag___('CSCG_form')
        self.standard_properties.___PRIVATE_add_tag___('form')


        self._CF_ = None
        self._BC_ = CSCG_Form_BC(self)

        self.___is_wrapped_in_ADF___ = False

    #------------------- fundamental --------------------------------
    @property
    def mesh(self):
        """Return the mesh."""
        return self._mesh_

    @property
    def space(self):
        """Return the basis function space."""
        return self._space_

    @property
    def p(self):
        """Return the degree of basis functions."""
        return self.space.p

    @property
    def k(self):
        return self._k_

    @property
    def ndim(self):
        """Return the dimensions."""
        return self.___ndim___


    @property
    def CF(self):
        """Continuous Form."""
        return self._CF_

    @CF.setter
    def CF(self, CF):
        self.___Pr_check_CF___(CF)
        self._CF_ = CF

    def ___Pr_check_CF___(self, CF):
        raise NotImplementedError()

    def ___Pr_check_BC_CF___(self, CF):
        raise NotImplementedError()





    @property
    def BC(self):
        return self._BC_

    # ------ must have properties ------------------------------------
    @property
    def cochain(self):
        raise NotImplementedError()

    # ------------ define information ----------------------------------
    @property
    def ___personal_parameters___(self):
        """
        Personal parameters are parameters that are additional to define_parameters.
        The full parameters
        are combination of both personal and define parameters.

        IMPORTANT:
        """

        # noinspection PyUnresolvedReferences
        RW_cochain_local = self.cochain.\
            ___PRIVATE_do_gather_to_master_and_make_them_region_wise_local_index_grouped___()
        # cochain_local will only be in master, in slaves, it is None.
        return {'region_wise_cochain_local': RW_cochain_local}

    @property
    def ___parameters___(self):
        """
        ___parameters___ will be called when you call standard_properties.parameters.

        :return:
        """
        parameters = dict()
        parameters.update(self.___define_parameters___)
        parameters.update(self.___personal_parameters___)
        return parameters