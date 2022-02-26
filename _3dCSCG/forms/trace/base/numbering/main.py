# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from screws.frozen import FrozenOnly
from importlib import import_module
from root.config import *





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
        path = '_3dCSCG.forms.trace.base.numbering.' + scheme_name
        name = '_3dCSCG_Trace_Form_Numbering_' + scheme_name
        self._numberer_ = getattr(import_module(path), name)(tf)
        self._parameters_ = parameters
        self._numbering_parameters_ = {'scheme_name': self._scheme_name_}
        self._numbering_parameters_.update(self._parameters_)
        self._DO_ = _3dCSCG_Trace_Numbering_DO(self)
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
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
    def num_of_dofs_in_this_core(self):
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

        RTE = boundaries.RANGE_trace_elements
        for bn in RTE:
            TE = RTE[bn]
            for te in TE:

                dofs = self.trace_element_wise[te]
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
        bd = cOmm.gather(LOCAL, root=mAster_rank)
        if rAnk == mAster_rank:
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
        return cOmm.bcast(BD, root=mAster_rank)

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
        """(dict) Return a dictionary of gathering vectors."""
        if self._trace_element_wise_ is None:
            self.___PRIVATE_do_numbering___()
        return self._trace_element_wise_


    @property
    def DO(self):
        return self._DO_







class _3dCSCG_Trace_Numbering_DO(FrozenOnly):
    def __init__(self, TN):
        self._numbering_ = TN
        self._find_ = _3dCSCG_Trace_Numbering_DO_FIND(self)
        self._freeze_self_()


    @property
    def FIND(self):
        return self._find_







class _3dCSCG_Trace_Numbering_DO_FIND(FrozenOnly):
    def __init__(self, DO):
        self._numbering_ = DO._numbering_
        self._0TraceLocalCache_ = dict()
        self._1TraceLocalCache_ = dict()
        self._2TraceLocalCache_ = dict()
        self._freeze_self_()

    def local_dofs_on_element_side(self, side_name):
        """

        :param side_name:
        :return:
        """
        numbering = self._numbering_

        if numbering._tf_.k == 0:
            return self.___PRIVATE_find_0Trace_local_dofs_on_element_side___(side_name)
        if numbering._tf_.k == 1:
            return self.___PRIVATE_find_1Trace_local_dofs_on_element_side___(side_name)
        if numbering._tf_.k == 2:
            return self.___PRIVATE_find_2Trace_local_dofs_on_element_side___(side_name)
        else:
            raise NotImplementedError(f"not coded for {numbering._tf_.k}-trace-form.")


    def ___PRIVATE_find_0Trace_local_dofs_on_element_side___(self, side_name):
        """

        :param side_name:
        :return:
        """
        if side_name in self._0TraceLocalCache_: return self._0TraceLocalCache_[side_name]

        nbc = self._numbering_._tf_.NUM_basis_components
        num_NS = nbc['N'][0]
        num_WE = nbc['W'][0]
        num_BF = nbc['B'][0]

        if side_name == 'N':
            self._0TraceLocalCache_['N'] = [i for i in range(num_NS)]
        elif side_name == 'S':
            self._0TraceLocalCache_['S'] = [i for i in range(num_NS, 2*num_NS)]
        elif side_name == 'W':
            self._0TraceLocalCache_['W'] = [i for i in range(2*num_NS, 2*num_NS+num_WE)]
        elif side_name == 'E':
            self._0TraceLocalCache_['E'] = [i for i in range(2*num_NS+num_WE, 2*num_NS+2*num_WE)]
        elif side_name == 'B':
            self._0TraceLocalCache_['B'] = [i for i in range(2*num_NS+2*num_WE, 2*num_NS+2*num_WE+num_BF)]
        elif side_name == 'F':
            self._0TraceLocalCache_['F'] = [i for i in range(2*num_NS+2*num_WE+num_BF, 2*num_NS+2*num_WE+2*num_BF)]
        else:
            raise Exception()

        return self._0TraceLocalCache_[side_name]

    def ___PRIVATE_find_1Trace_local_dofs_on_element_side___(self, side_name):
        """

        :param side_name:
        :return:
        """
        if side_name in self._1TraceLocalCache_: return self._1TraceLocalCache_[side_name]

        nbc = self._numbering_._tf_.NUM_basis_components
        num_NS = nbc['N'][0]
        num_WE = nbc['W'][0]
        num_BF = nbc['B'][0]

        if side_name == 'N':
            self._1TraceLocalCache_['N'] = [i for i in range(num_NS)]
        elif side_name == 'S':
            self._1TraceLocalCache_['S'] = [i for i in range(num_NS, 2*num_NS)]
        elif side_name == 'W':
            self._1TraceLocalCache_['W'] = [i for i in range(2*num_NS, 2*num_NS+num_WE)]
        elif side_name == 'E':
            self._1TraceLocalCache_['E'] = [i for i in range(2*num_NS+num_WE, 2*num_NS+2*num_WE)]
        elif side_name == 'B':
            self._1TraceLocalCache_['B'] = [i for i in range(2*num_NS+2*num_WE, 2*num_NS+2*num_WE+num_BF)]
        elif side_name == 'F':
            self._1TraceLocalCache_['F'] = [i for i in range(2*num_NS+2*num_WE+num_BF, 2*num_NS+2*num_WE+2*num_BF)]
        else:
            raise Exception()

        return self._1TraceLocalCache_[side_name]

    def ___PRIVATE_find_2Trace_local_dofs_on_element_side___(self, side_name):
        """

        :param side_name:
        :return:
        """
        if side_name in self._2TraceLocalCache_: return self._2TraceLocalCache_[side_name]

        nbc = self._numbering_._tf_.NUM_basis_components
        num_NS = nbc['N'][0]
        num_WE = nbc['W'][0]
        num_BF = nbc['B'][0]

        if side_name == 'N':
            self._2TraceLocalCache_['N'] = [i for i in range(num_NS)]
        elif side_name == 'S':
            self._2TraceLocalCache_['S'] = [i for i in range(num_NS, 2*num_NS)]
        elif side_name == 'W':
            self._2TraceLocalCache_['W'] = [i for i in range(2*num_NS, 2*num_NS+num_WE)]
        elif side_name == 'E':
            self._2TraceLocalCache_['E'] = [i for i in range(2*num_NS+num_WE, 2*num_NS+2*num_WE)]
        elif side_name == 'B':
            self._2TraceLocalCache_['B'] = [i for i in range(2*num_NS+2*num_WE, 2*num_NS+2*num_WE+num_BF)]
        elif side_name == 'F':
            self._2TraceLocalCache_['F'] = [i for i in range(2*num_NS+2*num_WE+num_BF, 2*num_NS+2*num_WE+2*num_BF)]
        else:
            raise Exception()

        return self._2TraceLocalCache_[side_name]