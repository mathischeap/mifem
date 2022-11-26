# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from components.freeze.main import FrozenOnly
from importlib import import_module
from objects.CSCG._3d.forms.standard.base.numbering.do.main import _3dCSCG_Standard_Form_Numbering_DO_


class _3dCSCG_Standard_Form_Numbering(FrozenOnly):
    """"""
    def __init__(self, sf, numbering_parameters):
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
        self._sf_ = sf
        self._scheme_name_ = scheme_name
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        path = base_path + scheme_name
        name = '_3dCSCG_Standard_Form_Numbering_' + scheme_name
        self._numberer_ = getattr(import_module(path), name)(sf)
        self._parameters_ = parameters
        self._numbering_parameters_ = {'scheme_name': self._scheme_name_}
        self._numbering_parameters_.update(self._parameters_)
        self._DO_ = _3dCSCG_Standard_Form_Numbering_DO_(self)
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        """"""
        self._local_ = None
        self._gathering_ = None
        self._boundary_dofs_ = None
        self._local_num_dofs_ = None
        self._extra_ = None
        self._localSideCache0_ = None
        self._localSideCache1_ = None
        self._localSideCache2_ = None

    def ___PRIVATE_do_numbering___(self):
        """"""
        if self._gathering_ is None or self._local_num_dofs_ is None or self._extra_ is None:
            self._gathering_, self._local_num_dofs_, self._extra_ = \
                getattr(self._numberer_, self._sf_.__class__.__name__)()

            assert len(self._gathering_) == len(self._sf_.mesh.elements), "GM length wrong."
            if SAFE_MODE:
                for i in self._sf_.mesh.elements: assert i in self._gathering_, "SAFETY test failed."
        else:
            pass

    @property
    def num_local_dofs(self):
        """"""
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
        """"""
        if self._gathering_ is None:
            self.___PRIVATE_do_numbering___()
        return self._gathering_

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
        mesh = self._sf_.mesh
        boundaries = mesh.boundaries
        BNS = boundaries.names
        self._boundary_dofs_ = dict()
        for bn in BNS: self._boundary_dofs_[bn] = list()
        if self._sf_.k == 3: return self._boundary_dofs_ # volume form has no boundary dofs!

        RES = boundaries.range_of_element_sides

        _ = self.gathering
        # this is important! Do so to avoid calling gathering in ``___PRIVATE_DO_find_dofs_on_element_side___``!

        for bn in RES:
            ES = RES[bn]
            for es in ES:
                e = int(es[:-1])
                s = es[-1]
                # dofs = self.___PRIVATE_DO_find_dofs_on_element_side___(e, s)
                dofs = self.do.find.dofs_on_element_side(e, s)
                # noinspection PyUnresolvedReferences
                self._boundary_dofs_[bn].extend(dofs)
            # noinspection PyUnresolvedReferences
            self._boundary_dofs_[bn] = list(set(self._boundary_dofs_[bn]))

        return self._boundary_dofs_

    @property
    def GLOBAL_boundary_dofs(self):
        """"""
        LOCAL = self.boundary_dofs
        if self._sf_.k == 3:
            return LOCAL
        else:
            bd = COMM.gather(LOCAL, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                mesh = self._sf_.mesh
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

                for name in names:
                    # noinspection PyUnresolvedReferences
                    BD[name] = list(set(BD[name]))
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
    def sharing_physical_locations(self):
        """Return a list of tuples which represent the dofs that locating at the same physical locations.

        For example, if it returns the following list:
            [(1, 45, 665, 128), (2, 88, 95), ...]
        Then we know that actually dofs #1, #45, #665, #128 exactly are at the same physical location.
        Similarly. Dofs #2, #88, #95 are at the same physical location. And so on.

        IMPORTANT: this is mainly for testing purpose. Its efficiency is very low. So for large simulations,
        calling this property will be extremely expensive.

        """
        if self._sf_.IS.hybrid is False: # for non-hybrid forms, no dofs are at the same location.
            return list()
        else:
            if self._sf_.k == 3: # a volume form, return empty list
                return list()
            else:
                kwargs = self._sf_.___define_parameters___['kwargs']
                kwargs['is_hybrid'] = False
                non_hybrid_form = self._sf_.__class__(self._sf_.mesh, self._sf_.space, **kwargs)
                NHF = non_hybrid_form.numbering.gathering
                NHF = COMM.gather(NHF, root=MASTER_RANK)

                GM = self.gathering
                GM = COMM.gather(GM, root=MASTER_RANK)
                if RANK == MASTER_RANK: # in the master rank.
                    DICT = dict()
                    for _, nhf in enumerate(NHF):
                        gm = GM[_]
                        for i in nhf: # we are going through all mesh elements.
                            GV = nhf[i]
                            gv = gm[i]
                            for j, dof in enumerate(GV):
                                h_dof = gv[j]
                                if dof not in DICT:
                                    DICT[dof] = (h_dof,)
                                else:
                                    DICT[dof] += (h_dof,)

                    LIST = list()
                    for i in DICT:
                        if len(DICT[i]) == 1:
                            pass
                        else:
                            LIST.append(DICT[i])

                    return LIST
                else:
                    return None

    @property
    def do(self):
        return self._DO_

    def ___PRIVATE_find_0Form_dofs_on_element_side___(self, element, side_name, GM):
        """

        :param element:
        :param side_name:
        :return:
        """
        assert self._sf_.num.basis == GM.GLOBAL_shape[1]
        if self._localSideCache0_ is None: self._localSideCache0_ = dict()
        if side_name not in self._localSideCache0_:
            if   side_name == 'N': indices = self.local[0][ 0, :, :]
            elif side_name == 'S': indices = self.local[0][-1, :, :]
            elif side_name == 'W': indices = self.local[0][ :, 0, :]
            elif side_name == 'E': indices = self.local[0][ :,-1, :]
            elif side_name == 'B': indices = self.local[0][ :, :, 0]
            elif side_name == 'F': indices = self.local[0][ :, :,-1]
            else: raise Exception()
            self._localSideCache0_[side_name] = indices.ravel('F')

        return GM[element].full_vector[self._localSideCache0_[side_name]]

    def ___PRIVATE_find_1Form_dofs_on_element_side___(self, element, side_name, GM):
        """

        :param element:
        :param side_name:
        :return:
        """
        assert self._sf_.num.basis == GM.GLOBAL_shape[1]
        if self._localSideCache1_ is None: self._localSideCache1_ = dict()
        if side_name not in self._localSideCache1_:
            if   side_name == 'N':
                indices1 = self.local[1][ 0, :, :]
                indices2 = self.local[2][ 0, :, :]
            elif side_name == 'S':
                indices1 = self.local[1][-1, :, :]
                indices2 = self.local[2][-1, :, :]
            elif side_name == 'W':
                indices1 = self.local[0][ :, 0, :]
                indices2 = self.local[2][ :, 0, :]
            elif side_name == 'E':
                indices1 = self.local[0][ :,-1, :]
                indices2 = self.local[2][ :,-1, :]
            elif side_name == 'B':
                indices1 = self.local[0][ :, :, 0]
                indices2 = self.local[1][ :, :, 0]
            elif side_name == 'F':
                indices1 = self.local[0][ :, :,-1]
                indices2 = self.local[1][ :, :,-1]
            else: raise Exception()
            self._localSideCache1_[side_name] = list()
            self._localSideCache1_[side_name].extend(indices1.ravel('F'))
            self._localSideCache1_[side_name].extend(indices2.ravel('F'))

        return GM[element].full_vector[self._localSideCache1_[side_name]]

    def ___PRIVATE_find_2Form_dofs_on_element_side___(self, element, side_name, GM):
        """

        :param element:
        :param side_name:
        :return:
        """
        assert self._sf_.num.basis == GM.GLOBAL_shape[1]
        if self._localSideCache2_ is None: self._localSideCache2_ = dict()
        if side_name not in self._localSideCache2_:
            if   side_name == 'N': indices = self.local[0][ 0, :, :]
            elif side_name == 'S': indices = self.local[0][-1, :, :]
            elif side_name == 'W': indices = self.local[1][ :, 0, :]
            elif side_name == 'E': indices = self.local[1][ :,-1, :]
            elif side_name == 'B': indices = self.local[2][ :, :, 0]
            elif side_name == 'F': indices = self.local[2][ :, :,-1]
            else: raise Exception()
            self._localSideCache2_[side_name] = indices.ravel('F')

        return GM[element].full_vector[self._localSideCache2_[side_name]]





if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\forms\standard\base\numbering\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy')([1, 1, 2])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    f = FC('2-f', is_hybrid=True)
    #
    bd = f.numbering.sharing_physical_locations

    print(bd)