# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from objects.base.form.base import FormBase
from objects.CSCG.base.forms.base.BC.main import CSCG_Form_BC

class CSCG_FORM_BASE(FormBase):
    """"""
    def __init_subclass__(cls, ndim=None):
        super().__init_subclass__()
        cls.___ndim___ = ndim

    def __init__(self, mesh, space, name):
        super(CSCG_FORM_BASE, self).__init__(mesh, space, name)
        self._k_ = None
        self.___define_parameters___ = None
        self.standard_properties.___PRIVATE_add_tag___('CSCG_form')

        self._BC_ = CSCG_Form_BC(self)

        self.___is_wrapped_in_ADF___ = False

    #------------------- fundamental --------------------------------
    @property
    def k(self):
        return self._k_

    @property
    def ndim(self):
        """Return the dimensions."""
        return self.___ndim___

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