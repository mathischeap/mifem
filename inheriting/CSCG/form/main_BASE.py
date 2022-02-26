# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from copy import deepcopy
from screws.frozen import FrozenClass
from inheriting.CSCG.form.time_wise import CSCG_Form_TimeWise
from inheriting.CSCG.form.time_wise import CSCG_Form_Func
from inheriting.CSCG.form.time_wise import CSCG_Form_BC_func

from tools.CSCG.partial_cochain import PartialCochain



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


    @property
    def TW(self):
        return self._TW_

    @property
    def func(self):
        return self._func_

    @property
    def BC(self):
        return self._BC_



    @property
    def cochain(self):
        raise NotImplementedError()

    @property
    def numbering(self):
        raise NotImplementedError()

    @property
    def visualize(self):
        raise NotImplementedError()



    @property
    def ___personal_parameters___(self):
        """
        Personal parameters are parameters that are additional to define_parameters. The full parameters
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
        RW_cochain_local = self.cochain.___PRIVATE_do_gather_to_master_and_make_them_region_wise_local_index_grouped___()
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
    def ndim(self):
        """Return the dimensions."""
        return self.___ndim___

    @property
    def k(self):
        return self._k_





class CSCG_Form_BC(CSCG_Form_BC_func):
    def __init__(self, f):
        super(CSCG_Form_BC, self).__init__(f)
        self._valid_boundaries_ = None # can not put it in RESET_cache method
        self.RESET_cache()
        self._freeze_self_()


    def RESET_cache(self):
        super(CSCG_Form_BC, self).RESET_cache()
        self._partial_cochain_ = None

    @property
    def valid_boundaries(self):
        return self._valid_boundaries_

    @valid_boundaries.setter
    def valid_boundaries(self, valid_boundaries):
        """This BC is valid on ``boundary_names``."""
        if isinstance(valid_boundaries, str):
            valid_boundaries = [valid_boundaries,]

        assert isinstance(valid_boundaries, (list, tuple)), \
            f"Please put boundary names in list or tuple"

        for i, bn in enumerate(valid_boundaries):
            assert bn in self._f_.mesh.boundaries.names, \
                f"boundary_names[{i}]: {bn} is not in mesh.boundaries.names: {self._f_.mesh.boundaries.names}"

        self._valid_boundaries_ = valid_boundaries
        self._partial_cochain_ = None


    @property
    def partial_cochain(self):
        """We will interpret the BC as a PartialCochain instance which then can
        further be interpreted as data structures that can be used by,
        for example, EWC sparse matrices.
        """
        if self._partial_cochain_ is None:
            pc = PartialCochain(self._f_)
            pc.include.boundaries(self.valid_boundaries)
            self._partial_cochain_ = pc
        return self._partial_cochain_