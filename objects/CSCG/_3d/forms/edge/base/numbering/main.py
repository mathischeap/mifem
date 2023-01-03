# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from components.freeze.main import FrozenOnly
from importlib import import_module
from objects.CSCG._3d.forms.edge.base.numbering.do import _3dCSCG_Edge_Numbering_DO


class _3dCSCG_Edge_Numbering(FrozenOnly):
    def __init__(self, ef, numbering_parameters):
        """"""
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
        self._ef_ = ef
        self._scheme_name_ = scheme_name
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        path = base_path + scheme_name
        name = '_3dCSCG_Edge_Numbering_' + scheme_name
        self._numberer_ = getattr(import_module(path), name)(ef)
        self._parameters_ = parameters
        self._numbering_parameters_ = {'scheme_name': self._scheme_name_}
        self._numbering_parameters_.update(self._parameters_)
        self._DO_ = _3dCSCG_Edge_Numbering_DO(self)
        self._gathering_ = None
        self._edge_element_wise_ = None
        self._local_num_dofs_ = None
        self._extra_ = None
        self._local_ = None
        self._boundary_dofs_ = None
        self._freeze_self_()

    @property
    def local(self):
        """The local numbering in mesh element."""
        if self._local_ is None:
            self._local_ = getattr(self._ef_.space.local_numbering, self._ef_.__class__.__name__)
        return self._local_

    def ___PRIVATE_do_numbering___(self):
        self._gathering_, self._edge_element_wise_, self._local_num_dofs_, self._extra_ = \
            getattr(self._numberer_, self._ef_.__class__.__name__)()

    @property
    def gathering(self):
        if self._gathering_ is None:
            self.___PRIVATE_do_numbering___()
        return self._gathering_

    @property
    def edge_element_wise(self):
        """(dict) Return a dictionary of gathering vectors."""
        if self._edge_element_wise_ is None:
            self.___PRIVATE_do_numbering___()
        return self._edge_element_wise_

    @property
    def num_local_dofs(self):
        if self._local_num_dofs_ is None:
            self.___PRIVATE_do_numbering___()
        return self._local_num_dofs_

    @property
    def extra(self):
        if self._extra_ is None:
            self.___PRIVATE_do_numbering___()
        return self._extra_

    @property
    def do(self):
        return self._DO_
