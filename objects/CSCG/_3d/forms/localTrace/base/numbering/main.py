# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/26/2022 2:27 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
from importlib import import_module
from root.config.main import RANK, MASTER_RANK, COMM, SAFE_MODE
from objects.CSCG._3d.forms.localTrace.base.numbering.do.main import _3dCSCG_LocalTrace_Numbering_Do

class _3dCSCG_LocalTrace_Numbering(FrozenOnly):
    """"""

    def __init__(self, ltf, numbering_parameters):
        """"""
        self._ltf_ = ltf

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

        self._scheme_name_ = scheme_name
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        path = base_path + scheme_name
        name = '_3dCSCG_LocalTrace_Form_Numbering_' + scheme_name
        self._numberer_ = getattr(import_module(path), name)(ltf)
        self._parameters_ = parameters
        self._numbering_parameters_ = {'scheme_name': self._scheme_name_}
        self._numbering_parameters_.update(self._parameters_)
        self._do_ = _3dCSCG_LocalTrace_Numbering_Do(ltf)
        self._local_gathering_ = None

        self._local_ = None

        self._gathering_ = None
        self._local_num_dofs_ = None
        self._extra_ = None

        self._boundary_dofs_ = None

        self._freeze_self_()

    @property
    def local_gathering(self):
        if self._ltf_.whether.hybrid:
            raise Exception(f"hybrid local-trace-form does not need local gathering matrix.")
        else:
            if self._local_gathering_ is None:
                self._local_gathering_ = getattr(
                    self._ltf_.space.local_gathering, self._ltf_.__class__.__name__
                )
            return self._local_gathering_

    @property
    def do(self):
        return self._do_

    @property
    def local(self):
        """The local numbering in mesh element."""
        if self._local_ is None:
            self._local_ = getattr(self._ltf_.space.local_numbering, self._ltf_.__class__.__name__)
        return self._local_

    def ___PRIVATE_do_numbering___(self):
        self._gathering_, self._local_num_dofs_, self._extra_ = \
            getattr(self._numberer_, self._ltf_.__class__.__name__)()

    @property
    def gathering(self):
        """The gathering matrix for this core."""
        if self._gathering_ is None:
            self.___PRIVATE_do_numbering___()
        return self._gathering_

    @property
    def num_local_dofs(self):
        """The amount of dofs in this core."""
        if self._local_num_dofs_ is None:
            self.___PRIVATE_do_numbering___()
        return self._local_num_dofs_

    @property
    def extra(self):
        """"""
        if self._extra_ is None:
            self.___PRIVATE_do_numbering___()
        return self._extra_


    @property
    def boundary_dofs(self):
        """local dofs on each boundary (exclude periodic boundaries).

        :return:
        """
        if self._boundary_dofs_ is not None:
            return self._boundary_dofs_
        mesh = self._ltf_.mesh
        boundaries = mesh.boundaries
        BNS = boundaries.names
        self._boundary_dofs_ = dict()
        for bn in BNS: self._boundary_dofs_[bn] = list()

        RES = boundaries.range_of_element_sides

        _ = self.gathering
        # this is important! Do so to avoid calling gathering during finding ``dofs_on_element_side``!

        for bn in RES:
            ES = RES[bn]
            for es in ES:
                e = int(es[:-1])
                s = es[-1]
                dofs = self.do.find.dofs_on_element_side(e, s)
                # noinspection PyUnresolvedReferences
                self._boundary_dofs_[bn].extend(dofs)
            # self._boundary_dofs_[bn] = list(set(self._boundary_dofs_[bn])) # no need to do this as it is hybrid.

        return self._boundary_dofs_

    @property
    def global_boundary_dofs(self):
        """Return all dofs on boundary in a dict whose keys are boundary names and values are dof
        numbers in all cores."""
        LOCAL = self.boundary_dofs

        bd = COMM.gather(LOCAL, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            mesh = self._ltf_.mesh
            BD = dict()
            names = mesh.boundaries.names
            for name in names:
                BD[name] = list()
            for bdi in bd:
                for name in names:
                    # noinspection PyUnresolvedReferences
                    A = BD[name]
                    if SAFE_MODE: assert isinstance(A, list), "SAFETY check failed."
                    A.extend(bdi[name])
                    BD[name] = A

            # for name in names:                   # no need as it is hybrid.
            #     BD[name] = list(set(BD[name]))
        else:
            BD = None
        return COMM.bcast(BD, root=MASTER_RANK)

    @property
    def global_boundary_dofs_ravel(self):
        """Regardless to boundary names; collection of all boundary dofs in one set."""
        GLOBAL_boundary_dofs = self.global_boundary_dofs
        RAVEL = set()
        for bn in GLOBAL_boundary_dofs:
            RAVEL.update(GLOBAL_boundary_dofs[bn])
        return list(RAVEL)

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
