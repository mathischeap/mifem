# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from screws.freeze.main import FrozenOnly
from importlib import import_module






class _3dCSCG_Tr_Numbering(FrozenOnly):
    def __init__(self, Tr, numbering_parameters):
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

        self._Tr_ = Tr
        self._scheme_name_ = scheme_name
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        path = base_path + scheme_name
        name = '_3dCSCG_TrForm_Numbering_' + scheme_name
        self._numberer_ = getattr(import_module(path), name)(Tr)
        self._parameters_ = parameters
        self._numbering_parameters_ = {'scheme_name': self._scheme_name_}
        self._numbering_parameters_.update(self._parameters_)

        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._local_ = None
        self._local_gathering_ = None
        self._gathering_ = None
        self._trace_element_wise_ = None
        self._boundary_dofs_ = None
        self._local_num_dofs_ = None
        self._extra_ = None

    @property
    def local(self):
        """The local numbering in mesh element."""
        if self._local_ is None:
            self._local_, self._local_gathering_ = getattr(
                self._Tr_.space.local_numbering, self._Tr_.__class__.__name__
            )
        return self._local_

    @property
    def local_gathering(self):
        """The local numbering in mesh element."""
        if self._local_gathering_ is None:
            self._local_, self._local_gathering_ = getattr(
                self._Tr_.space.local_numbering, self._Tr_.__class__.__name__
            )
        return self._local_gathering_





    def ___PRIVATE_do_numbering___(self):
        self._gathering_, self._trace_element_wise_, self._local_num_dofs_, self._extra_ = \
            getattr(self._numberer_, self._Tr_.__class__.__name__)()

    @property
    def gathering(self):
        if self._gathering_ is None:
            self.___PRIVATE_do_numbering___()
        return self._gathering_

    @property
    def trace_element_wise(self):
        """(dict) Return a dictionary of trace-element-wise gathering vectors."""
        if self._trace_element_wise_ is None:
            self.___PRIVATE_do_numbering___()
        return self._trace_element_wise_

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
    def boundary_dofs(self):
        """Return a dict whose keys are all (in all cores) mesh boundary names,
        and values are the local dofs on the boundary indicated by the key.
        """
        if self._boundary_dofs_ is not None:
            return self._boundary_dofs_
        mesh = self._Tr_.mesh
        boundaries = mesh.boundaries
        BNS = boundaries.names
        self._boundary_dofs_ = dict()
        for bn in BNS: self._boundary_dofs_[bn] = list()

        RTE = boundaries.range_of_trace_elements
        for bn in RTE:
            TE = RTE[bn]
            for te in TE:

                dofs = self.trace_element_wise[te]
                # noinspection PyUnresolvedReferences
                self._boundary_dofs_[bn].extend(dofs)

            # noinspection PyUnresolvedReferences
            self._boundary_dofs_[bn] = list(set(self._boundary_dofs_[bn]))

        return self._boundary_dofs_