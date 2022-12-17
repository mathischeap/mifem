# -*- coding: utf-8 -*-
import sys
if './' not in sys.path:
    sys.path.append('./')

from root.config.main import *
from importlib import import_module
from components.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.domain.inputs.allocator import DomainInputFinder
from objects.CSCG._2d.mesh.domain.main import _2dCSCG_Domain
from objects.CSCG._2d.mesh.main import _2dCSCG_Mesh

from objects.CSCG._2d.spaces.allocator import _2dCSCG_SpaceAllocator

from objects.CSCG._2d.forms.allocator import _2dCSCG_FormsAllocator
from objects.CSCG._2d.fields.allocator import _2dCSCG_FieldsAllocator

from objects.CSCG._2d.exactSolutions.allocator import _2dCSCG_ExactSolutionAllocator

from copy import deepcopy
from components.miscellaneous.timer import MyTimer


class MeshGenerator(FrozenOnly):
    def __init__(self, ID, **kwargs):
        """Remember, **kwargs are parameters to customize the domain.
        The rule is: they can not change the topology of the regions!
        """
        COMM.barrier()  # for safety reason
        di = DomainInputFinder(ID)(**kwargs)
        self._domain_ = _2dCSCG_Domain(di)
        self._ID_ = ID
        self._kwargs_ = kwargs
        self._freeze_self_()

    def __call__(self, element_layout, EDM=None, show_info=False):
        if show_info and RANK == MASTER_RANK:
            print(f"---[2dCSCG]-[MESH]-{MyTimer.current_time()}-----")
            print(f"   <domain ID>: {self._ID_}")
            str_dp = str(self._kwargs_)
            if len(str_dp) > 40: str_dp = str_dp[:40] + '...'
            print( "   <domain_parameters>: {}".format(str_dp))
            print(f"   <EDM>: {EDM}", flush=True)

        COMM.barrier()  # for safety reason
        mesh = _2dCSCG_Mesh(self._domain_, element_layout, EDM=EDM)
        mp = dict()
        mp['type'] = '_2dCSCG_Mesh'
        mp['ID'] = self._ID_
        mp['domain_parameters'] = deepcopy(self._kwargs_)
        mp['element_layout'] = element_layout
        mp['EDM'] = EDM
        mesh.___define_parameters___ = mp
        dp = dict()
        dp['ID'] = self._ID_
        for key in self._kwargs_:
            dp[key] = self._kwargs_[key]
        mesh.domain.___define_parameters___ = dp

        if show_info and RANK == MASTER_RANK:
            str_element_layout = str(element_layout)
            if len(str_element_layout) < 40:
                print( "   <element_layout input>: {}".format(str_element_layout))
            else:
                print( "   <element_layout input>: {}...".format(str_element_layout[:40]))
            for rn in mesh.elements.layout:
                print(f"   <element_layout>: {rn} {mesh.elements.layout[rn]}")
            print(f"   <total elements>: {mesh.elements.global_num}", flush=True)

        COMM.barrier()  # for safety reason

        return mesh

    @classmethod
    def ___coded_meshes___(cls):
        return DomainInputFinder.___defined_DI___()

    @classmethod
    def ___domain_input_statistic___(cls):
        Statistic = dict()
        for diID in DomainInputFinder.___defined_DI___():
            cls_name = DomainInputFinder.___defined_DI___()[diID]
            cls_path = DomainInputFinder.___DI_path___()[diID]
            diClass = getattr(import_module(cls_path), cls_name)
            try:
                Statistic[diID] = diClass.statistic
            except NotImplementedError:
                pass
        return Statistic

    @classmethod
    def ___domain_input_random_parameters___(cls):
        random_parameters = dict()
        for diID in DomainInputFinder.___defined_DI___():
            cls_name = DomainInputFinder.___defined_DI___()[diID]
            cls_path = DomainInputFinder.___DI_path___()[diID]
            diClass = getattr(import_module(cls_path), cls_name)
            try:
                random_parameters[diID] = diClass.random_parameters
            except NotImplementedError:
                pass
        return random_parameters


class SpaceInvoker(FrozenOnly):
    def __init__(self, ID):
        COMM.barrier()  # for safety reason
        assert ID in _2dCSCG_SpaceAllocator.___space_name___(), \
            " <SpaceInvoker> : space <{}> is not coded yet.".format(ID)
        self._ID_ = ID
        cls_name = _2dCSCG_SpaceAllocator.___space_name___()[ID]
        cls_path = _2dCSCG_SpaceAllocator.___space_path___()[ID]
        self._space_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, inputs, ndim=None, show_info=False):
        if show_info and RANK == MASTER_RANK:
            print(f"---[2dCSCG]-[SPACE]-{MyTimer.current_time()}-----")
            print(f"   <space ID>: {self._ID_}")
            print(f"   <space inputs>: {inputs}")

        COMM.barrier()  # for safety reason
        if ndim is not None: assert ndim == 2
        S = self._space_(inputs, ndim)
        sp = dict()
        sp['type'] = '_2dCSCG_Space'
        sp['ID'] = self._ID_
        sp['inputs'] = inputs
        sp['ndim'] = ndim
        S.___define_parameters___ = sp
        COMM.barrier()  # for safety reason
        return S


class FormCaller(FrozenOnly):
    """Generate a form."""
    def __init__(self, mesh, space):
        COMM.barrier()  # for safety reason
        self._mesh_ = mesh
        self._space_ = space

        PATH = dict()
        NAME = dict()
        for A in (_2dCSCG_FormsAllocator, _2dCSCG_FieldsAllocator):
            PATH.update(A.___form_path___())
            NAME.update(A.___form_name___())
        self._PATH_ = PATH
        self._NAME_ = NAME

        self._freeze_self_()

    def __call__(self, ID, *args, **kwargs):
        COMM.barrier()  # for safety reason

        cls_body = getattr(import_module(self._PATH_[ID]), self._NAME_[ID])
        if ID in ('scalar', 'vector'):
            FM = cls_body(self._mesh_, *args, **kwargs)
            # continuous forms has no parameters, not a FrozenClass. SO we CAN NOT save them.

        else:
            FM = cls_body(self._mesh_, self._space_, **kwargs)
            fp = dict()
            fp['type'] = '_2dCSCG_Form'
            fp['ID'] = ID
            fp['mesh_parameters'] = deepcopy(self._mesh_.standard_properties.parameters)
            fp['space_parameters'] = deepcopy(self._space_.standard_properties.parameters)
            fp['kwargs'] = kwargs
            FM.___define_parameters___ = fp
        COMM.barrier()  # for safety reason
        return FM

    @property
    def mesh(self):
        return self._mesh_

    @property
    def space(self):
        return self._space_


class ExactSolutionSelector(FrozenOnly):
    """We select an exact solution object with this class."""
    def __init__(self, mesh):
        COMM.barrier()  # for safety reason
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, ID, show_info=False, **kwargs):
        if show_info and RANK == MASTER_RANK:
            print(f"---[2dCSCG]-[Exact Solution]-{MyTimer.current_time()}-----")
            print(f"   <ES ID>: {ID}")
            print(f"   <ES kwargs>: {kwargs}", flush=True)

        COMM.barrier()  # for safety reason
        assert ID in _2dCSCG_ExactSolutionAllocator.___exact_solution_name___(), \
            f"Exact solution ID={ID} not found."
        className = _2dCSCG_ExactSolutionAllocator.___exact_solution_name___()[ID]
        classPath = _2dCSCG_ExactSolutionAllocator.___exact_solution_path___()[ID]

        ES = getattr(import_module(classPath), className)(self._mesh_, **kwargs)

        return ES

    @classmethod
    def listing(cls, printing=True, returning=False):
        """For an allocator class, this lists all the possibilities ONLY in the master core."""
        if RANK != MASTER_RANK: return
        included_allocators = [
            _2dCSCG_ExactSolutionAllocator,
        ]
        listing = ''
        for alc in included_allocators:
            listing += alc.listing(printing=False, returning=True)

        if printing:
            print(listing)
        else:
            pass
        if returning:
            return listing
        else:
            pass


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/master.py
    # import numpy as np

    # mesh = MeshGenerator('cic',)([14,14])
    # mesh = MeshGenerator('rectangle', p_UL=(-1,-1),region_layout=(3,5))([5,5], show_info=True)
    # mesh = MeshGenerator('rectangle_periodic', p_UL=(-1,-1),region_layout=(3,5))([5,5], show_info=True)
    # mesh = MeshGenerator('rectangle')([1,1], show_info=True)

    ExactSolutionSelector.listing()
