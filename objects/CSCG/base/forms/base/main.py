# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from copy import deepcopy
from screws.freeze.main import FrozenClass

from objects.CSCG.base.forms.base.time_wise.main import CSCG_Form_TimeWise
from objects.CSCG.base.forms.base.func import CSCG_Form_Func
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

        self._TW_ = CSCG_Form_TimeWise(self) # TW has 1) func
        self._func_ = CSCG_Form_Func(self)
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

    # --------------- base properties ---------------------------------
    @property
    def TW(self):
        return self._TW_

    @property
    def func(self):
        return self._func_

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
        if self.TW.func.body is None:
            ES_parameters = None
        else:
            if self.TW.func._ES_ is None:
                ES_parameters = None
            else:
                ES_parameters = \
                    [deepcopy(self.TW.func._ES_.standard_properties.parameters),
                     self.TW.func._ES_variable_name_]
        # noinspection PyUnresolvedReferences
        RW_cochain_local = self.cochain.\
            ___PRIVATE_do_gather_to_master_and_make_them_region_wise_local_index_grouped___()
        # cochain_local will only be in master, in slaves, it is None.
        return {'TW_current_time': self.TW.current_time, # to save a object, it must have current time.
                'TW_func_ES_parameters': ES_parameters,
                'region_wise_cochain_local': RW_cochain_local}

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