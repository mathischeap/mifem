# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from SCREWS.frozen import FrozenOnly
from importlib import import_module
from _2dCSCG.form.standard.numbering.visualize import _2dCSCG_Numbering_Visualize



class _2dCSCG_Standard_Form_Numbering(FrozenOnly):
    def __init__(self, sf, numbering_parameters):
        # ... parse number and numbering parameters ...
        if isinstance(numbering_parameters, str):
            scheme_name = numbering_parameters
            parameters = dict()
        elif isinstance(numbering_parameters, dict):
            scheme_name = numbering_parameters['scheme_name']
            parameters = dict()
            for key in numbering_parameters:
                if key != 'scheme_name':
                    parameters[key] = numbering_parameters[key]
        else:
            raise NotImplementedError()
        # ...
        self._sf_ = sf
        self._scheme_name_ = scheme_name
        path = '_2dCSCG.form.standard.numbering.' + scheme_name
        name = '_2dCSCG_Standard_Form_Numbering_' + scheme_name
        self._numberer_ = getattr(import_module(path), name)(sf)
        self._parameters_ = parameters
        self._numbering_parameters_ = {'scheme_name': self._scheme_name_}
        self._numbering_parameters_.update(self._parameters_)
        self._visualize_ = _2dCSCG_Numbering_Visualize(sf)
        self.DO_reset_cache()
        self._freeze_self_()

    def DO_reset_cache(self):
        self._local_ = None
        self._gathering_ = None
        self._boundary_dofs_ = None
        self._local_num_dofs_ = None
        self._extra_ = None

    def ___PRIVATE_do_numbering___(self):
        self._gathering_, self._local_num_dofs_, self._extra_ = \
            getattr(self._numberer_, self._sf_.__class__.__name__)()

    @property
    def num_of_dofs_in_this_core(self):
        if self._local_num_dofs_ is None:
            self.___PRIVATE_do_numbering___()
        return self._local_num_dofs_

    @property
    def local(self):
        """The local numbering in mesh element."""
        if self._local_ is None:
            self._local_ = getattr(self._sf_.space.local_numbering, self._sf_.__class__.__name__)
        return self._local_

    @property
    def gathering(self):
        if self._gathering_ is None:
            self.___PRIVATE_do_numbering___()
        return self._gathering_

    @property
    def boundary_dofs(self):
        raise NotImplementedError(f"Boundary dofs are not coded in numbering scheme yet.")

    @property
    def extra(self):
        if self._extra_ is None:
            self.___PRIVATE_do_numbering___()
        return self._extra_

    @property
    def visualize(self):
        return self._visualize_