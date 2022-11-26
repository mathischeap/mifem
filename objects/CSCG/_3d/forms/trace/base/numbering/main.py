# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from components.freeze.main import FrozenOnly
from importlib import import_module
from root.config.main import *


from objects.CSCG._3d.forms.trace.base.numbering.do.main import _3dCSCG_Trace_Numbering_DO



class _3dCSCG_Trace_Numbering(FrozenOnly):
    def __init__(self, tf, numbering_parameters):
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
        self._tf_ = tf
        self._scheme_name_ = scheme_name
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        path = base_path + scheme_name
        name = '_3dCSCG_Trace_Form_Numbering_' + scheme_name
        self._numberer_ = getattr(import_module(path), name)(tf)
        self._parameters_ = parameters
        self._numbering_parameters_ = {'scheme_name': self._scheme_name_}
        self._numbering_parameters_.update(self._parameters_)
        self._DO_ = _3dCSCG_Trace_Numbering_DO(self)
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._local_ = None
        self._gathering_ = None
        self._trace_element_wise_ = None
        self._boundary_dofs_ = None
        self._local_num_dofs_ = None
        self._extra_ = None

    def ___PRIVATE_do_numbering___(self):
        self._gathering_, self._trace_element_wise_, self._local_num_dofs_, self._extra_ = \
            getattr(self._numberer_, self._tf_.__class__.__name__)()

    @property
    def num_local_dofs(self):
        if self._local_num_dofs_ is None:
            self.___PRIVATE_do_numbering___()
        return self._local_num_dofs_

    @property
    def local(self):
        """The local numbering in mesh element."""
        if self._local_ is None:
            self._local_ = getattr(self._tf_.space.local_numbering, self._tf_.__class__.__name__)
        return self._local_

    @property
    def gathering(self):
        if self._gathering_ is None:
            self.___PRIVATE_do_numbering___()
        return self._gathering_

    @property
    def boundary_dofs(self):
        """Return a dict whose keys are all (in all cores) mesh boundary names,
        and values are the local dofs on the boundary indicated by the key.
        """
        if self._boundary_dofs_ is not None:
            return self._boundary_dofs_
        mesh = self._tf_.mesh
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

    @property
    def GLOBAL_boundary_dofs(self):
        """Return in all cores. Return a dict whose keys are mesh
        boundary names and values are all (in all cores) dofs on the boundary indicated by the key.
        """
        LOCAL = self.boundary_dofs
        bd = COMM.gather(LOCAL, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            mesh = self._tf_.mesh
            BD = dict()
            names = mesh.boundaries.names
            for name in names:
                BD[name] = list()
            for bdi in bd:
                for name in names:
                    # noinspection PyUnresolvedReferences
                    BD[name].extend(bdi[name])
            for name in names:
                # noinspection PyUnresolvedReferences
                BD[name] = set(BD[name])
        else:
            BD = None
        return COMM.bcast(BD, root=MASTER_RANK)

    @property
    def GLOBAL_boundary_dofs_ravel(self):
        """Regardless to boundary names; collection of all boundary dofs in one set."""
        GLOBAL_boundary_dofs = self.GLOBAL_boundary_dofs
        RAVEL = set()
        for bn in GLOBAL_boundary_dofs:
            RAVEL.update(GLOBAL_boundary_dofs[bn])
        return list(RAVEL)

    @property
    def extra(self):
        if self._extra_ is None:
            self.___PRIVATE_do_numbering___()
        return self._extra_

    @property
    def trace_element_wise(self):
        """(dict) Return a dictionary of trace-element-wise gathering vectors."""
        if self._trace_element_wise_ is None:
            self.___PRIVATE_do_numbering___()
        return self._trace_element_wise_


    @property
    def do(self):
        return self._DO_

