# -*- coding: utf-8 -*-
"""INTRO

@author: Yi Zhang. Created on Thu Jan 23 00:00:12 2020
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands
"""
import sys
if './' not in sys.path: sys.path.append('./')

from importlib import import_module
from screws.freeze.main import FrozenOnly
from screws.miscellaneous.timer import MyTimer
from objects.CSCG._3d.mesh.domain.inputs.allocator import DomainInputAllocator
from objects.CSCG._3d.mesh.domain.main import _3dCSCG_Domain
from objects.CSCG._3d.mesh.main import _3dCSCG_Mesh
from objects.CSCG._3d.spaces.allocator import _3dCSCG_SpaceAllocator
from objects.CSCG._3d.ADF.allocator import _3dCSCG_ADF_Allocator
from objects.CSCG._3d.forms.allocator import _3dCSCG_SF_Allocator
from objects.CSCG._3d.fields.allocator import _3dCSCG_Field_Allocator
from objects.CSCG._3d.exactSolutions.allocator import _3dCSCG_ExactSolution_Allocator
from copy import deepcopy
from root.config.main import RANK, MASTER_RANK, COMM


class MeshGenerator(FrozenOnly):
    def __init__(self, ID, **kwargs):
        """Remember, **kwargs are parameters to customize the domain.
        The rule is: they can not change the topology of the regions!
        """
        di = DomainInputAllocator(ID)(**kwargs)
        self._domain_ = _3dCSCG_Domain(di)
        self._ID_ = ID
        self._kwargs_ = kwargs
        self._freeze_self_()

    def __call__(self, element_layout, EDM=None, show_info=False):
        """
        :param element_layout:
        :param EDM: ``Element Distribution Method``. When
            EDM = 'debug', we use a naive method
            EDM = None, we find a proper one.
            EDM = a string, we use the specific one.
        :param show_info:
        :return:
        """
        if show_info and RANK == MASTER_RANK:
            print(f"---[3dCSCG]-[MESH]-{MyTimer.current_time()}-----")
            print(f"   <domain ID>: {self._ID_}")
            str_dp = str(self._kwargs_)
            if len(str_dp) > 40: str_dp = str_dp[:40] + '...'
            print( "   <domain_parameters>: {}".format(str_dp))
            print(f"   <EDM>: {EDM}", flush=True)

        COMM.barrier()  # for safety reason
        mesh = _3dCSCG_Mesh(self._domain_, element_layout, EDM=EDM)
        mp = dict()
        mp['type'] = '_3dCSCG_Mesh'
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
            print(f"   <total elements>: {mesh.elements.GLOBAL_num}", flush=True)
        COMM.barrier()  # for safety reason

        return mesh

    @classmethod
    def ___coded_meshes___(cls):
        return DomainInputAllocator.___defined_DI___()

    @classmethod
    def ___domain_input_statistic___(cls):
        Statistic = dict()
        for diID in DomainInputAllocator.___defined_DI___():
            cls_name = DomainInputAllocator.___defined_DI___()[diID]
            cls_path = DomainInputAllocator.___DI_path___()[diID]
            diClass = getattr(import_module(cls_path), cls_name)
            try:
                Statistic[diID] = diClass.statistic
            except NotImplementedError:
                pass
        return Statistic

    @classmethod
    def ___domain_input_random_parameters___(cls):
        random_parameters = dict()
        for diID in DomainInputAllocator.___defined_DI___():
            cls_name = DomainInputAllocator.___defined_DI___()[diID]
            cls_path = DomainInputAllocator.___DI_path___()[diID]
            diClass = getattr(import_module(cls_path), cls_name)
            try:
                random_parameters[diID] = diClass.random_parameters
            except NotImplementedError:
                pass
        return random_parameters


class SpaceInvoker(FrozenOnly):
    def __init__(self, ID):
        COMM.barrier()  # for safety reason
        assert ID in _3dCSCG_SpaceAllocator.___space_name___(), \
            " <SpaceInvoker> : space <{}> is not coded yet.".format(ID)
        self._ID_ = ID
        cls_name = _3dCSCG_SpaceAllocator.___space_name___()[ID]
        cls_path = _3dCSCG_SpaceAllocator.___space_path___()[ID]
        self._space_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, inputs, ndim=None, show_info=False):
        if show_info and RANK == MASTER_RANK:
            print(f"---[3dCSCG]-[SPACE]-{MyTimer.current_time()}-----")
            print(f"   <space ID>: {self._ID_}")
            print(f"   <space inputs>: {inputs}")
            print(f"   <space ndim>: {ndim}", flush=True)

        COMM.barrier()  # for safety reason
        S = self._space_(inputs, ndim)
        sp = dict()
        sp['type'] = '_3dCSCG_Space'
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
        for A in (_3dCSCG_ADF_Allocator, _3dCSCG_SF_Allocator, _3dCSCG_Field_Allocator):
            PATH.update(A.___forms_path___())
            NAME.update(A.___forms_name___())
        self._PATH_ = PATH
        self._NAME_ = NAME
        self._freeze_self_()

    def __call__(self, ID, *args, **kwargs):
        COMM.barrier()  # for safety reason

        cls_body = getattr(import_module(self._PATH_[ID]), self._NAME_[ID])

        if ID in ('scalar', 'vector', 'tensor'):

            FM = cls_body(self._mesh_, *args, **kwargs)
            # We CANNOT (do NOT) save continuous field instances.

        else:
            fp = dict()
            fp['ID'] = ID
            fp['mesh_parameters'] = deepcopy(self._mesh_.standard_properties.parameters)
            fp['space_parameters'] = deepcopy(self._space_.standard_properties.parameters)
            fp['kwargs'] = deepcopy(kwargs)

            if ID in ('0-adf', '1-adf', '2-adf', '3-adf',  # algebraic dual (standard) forms
                      '0-adt', '1-adt', '2-adt',           # algebraic dual trace forms
                      ):
                fp['type'] = '_3dCSCG_ADF'
                # ---------------- make a dual from a prime ------------------------------------
                if len(args) == 1: # if so, we get a prime form, we make dual form from it.
                    assert kwargs == dict(), \
                        f"when make algebraic dual standard form from prime form, " \
                        f"kwargs must be empty, now it is={kwargs}."
                    prime = args[0]
                    pcn = prime.__class__.__name__

                    if ID == '0-adf':
                        assert pcn == '_3dCSCG_0Form'
                    elif ID == '1-adf':
                        assert pcn == '_3dCSCG_1Form'
                    elif ID == '2-adf':
                        assert pcn == '_3dCSCG_2Form'
                    elif ID == '3-adf':
                        assert pcn == '_3dCSCG_3Form'
                    elif ID == '0-adt':
                        assert pcn == '_3dCSCG_0Trace'
                    elif ID == '1-adt':
                        assert pcn == '_3dCSCG_1Trace'
                    elif ID == '2-adt':
                        assert pcn == '_3dCSCG_2Trace'
                    else:
                        raise Exception(f"ID={ID} do not accept a single prime form instance as input.")

                    assert prime.IS.hybrid, "prime must be hybrid."
                    assert prime.mesh is self._mesh_ and prime.space is self._space_, \
                        "mesh, space do not match." # not just ==, but is!

                    # if we made a dual from a prime, we still need to feed `orientation` and `name` to the dual class.
                    kwargs = dict()
                    kwargs['orientation'] = prime.orientation
                    kwargs['name'] = 'AD_' + prime.standard_properties.name

                #--------- make the prime from the args and kwargs ---------------------------------
                else: # we make the prime form from the args.
                    p_kwargs = dict()
                    for kw in kwargs:
                        if kw == 'name':
                            p_kwargs[kw] = 'PRIME_' +  kwargs[kw]
                        else:
                            p_kwargs[kw] = kwargs[kw]

                    if ID in ('0-adf', '1-adf', '2-adf', '3-adf'):
                        prime_class_ID = ID.split('-')[0] + '-f'
                        prime = self(prime_class_ID, *args, **p_kwargs, is_hybrid=True)

                    elif ID in ('0-adt', '1-adt', '2-adt'): # note that trace forms must be hybrid.
                        prime_class_ID = ID.split('-')[0] + '-t'
                        prime = self(prime_class_ID, *args, **p_kwargs)

                    else:
                        raise Exception()

                # ---- remove `numbering_parameters` for dual class initialize -------------------
                KWARGS = dict()
                for key in kwargs: # we remove `numbering_parameters` form the kwargs because we don't need it
                                   # form the initialization of the dual form.
                    if key == 'numbering_parameters':
                        pass
                    else:
                        KWARGS[key] = kwargs[key]
                # == ABOVE =========================================================================

                FM = cls_body(prime, self._mesh_, self._space_, **KWARGS)

            else: # all other forms, like the standard forms, trace forms, edge forms, node forms and so on
                assert len(args) == 0, "all these forms do not take args, only take kwargs."
                fp['type'] = '_3dCSCG_Form'
                FM = cls_body(self._mesh_, self._space_, **kwargs)

            #====================================================================================
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
        """"""
        if show_info and RANK == MASTER_RANK:
            print(f"---[3dCSCG]-[Exact Solution]-{MyTimer.current_time()}-----")
            print(f"   <ES ID>: {ID}")
            print(f"   <ES kwargs>: {kwargs}", flush=True)

        COMM.barrier()  # for safety reason
        assert ID in _3dCSCG_ExactSolution_Allocator.___exact_solution_name___(), \
            f"Exact solution ID={ID} not found."
        className = _3dCSCG_ExactSolution_Allocator.___exact_solution_name___()[ID]
        classPath = _3dCSCG_ExactSolution_Allocator.___exact_solution_path___()[ID]

        ES = getattr(import_module(classPath), className)(self._mesh_, **kwargs)

        return ES





if __name__ == "__main__":
    # mpiexec -n 8 python objects\CSCG\_3d\master.py

    # mesh = MeshGenerator('cuboid', region_layout=(3,3,3))([2,3,4], show_info=True)
    # mesh = MeshGenerator('bridge_arch_cracked')([2,2,2], EDM='SWV0', show_info=True)
    mesh = MeshGenerator('crazy')([3, 3, 3], EDM='chaotic', show_info=True)
    # mesh = MeshGenerator('crazy_periodic')([3, 3, 3], EDM='chaotic', show_info=True)
    # mesh = MeshGenerator('cuboid_periodic', region_layout=(1,2,2))([2,3,4], show_info=True)

    # es = ExactSolutionSelector(mesh)('TISE:sincos1', m=1e-68, E=1)

    es = ExactSolutionSelector(mesh)('MHD:sincos1')
