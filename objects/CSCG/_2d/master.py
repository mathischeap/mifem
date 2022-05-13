# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from objects.CSCG._2d.exact_solutions.main import ExactSolution
from importlib import import_module
from screws.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.domain.inputs.allocator import DomainInputFinder
from objects.CSCG._2d.mesh.domain.main import _2dCSCG_Domain
from objects.CSCG._2d.mesh.main import _2dCSCG_Mesh

from objects.CSCG._2d.spaces.allocator import _2dCSCG_SpaceAllocator

from objects.CSCG._2d.forms.allocator import _2dCSCG_FormsAllocator
from objects.CSCG._2d.fields.allocator import _2dCSCG_FieldsAllocator

from objects.CSCG._2d.exact_solutions.status.allocator import _2dCSCG_ExactSolutionAllocator

from copy import deepcopy
from screws.miscellaneous.timer import MyTimer



class MeshGenerator(FrozenOnly):
    def __init__(self, ID, **kwargs):
        """Remember, **kwargs are parameters to customize the domain.
        The rule is: they can not change the topology of the regions!
        """
        cOmm.barrier()  # for safety reason
        di = DomainInputFinder(ID)(**kwargs)
        self._domain_ = _2dCSCG_Domain(di)
        self._ID_ = ID
        self._kwargs_ = kwargs
        self._freeze_self_()


    def __call__(self, element_layout, EDM=None, show_info=False):
        if show_info and rAnk == mAster_rank:
            print(f"---[2dCSCG]-[MESH]-{MyTimer.current_time()}-----")
            print(f"   <domain ID>: {self._ID_}")
            str_dp = str(self._kwargs_)
            if len(str_dp) > 40: str_dp = str_dp[:40] + '...'
            print( "   <domain_parameters>: {}".format(str_dp))
            print(f"   <EDM>: {EDM}", flush=True)

        cOmm.barrier()  # for safety reason
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

        if show_info and rAnk == mAster_rank:
            str_element_layout = str(element_layout)
            if len(str_element_layout) < 40:
                print( "   <element_layout input>: {}".format(str_element_layout))
            else:
                print( "   <element_layout input>: {}...".format(str_element_layout[:40]))
            for rn in mesh.elements.layout:
                print(f"   <element_layout>: {rn} {mesh.elements.layout[rn]}")
            print(f"   <total elements>: {mesh.elements.GLOBAL_num}", flush=True)

        cOmm.barrier()  # for safety reason

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
            Statistic[diID] = diClass.statistic
        return Statistic

    @classmethod
    def ___domain_input_random_parameters___(cls):
        random_parameters = dict()
        for diID in DomainInputFinder.___defined_DI___():
            cls_name = DomainInputFinder.___defined_DI___()[diID]
            cls_path = DomainInputFinder.___DI_path___()[diID]
            diClass = getattr(import_module(cls_path), cls_name)
            random_parameters[diID] = diClass.random_parameters
        return random_parameters





class SpaceInvoker(FrozenOnly):
    def __init__(self, ID):
        cOmm.barrier()  # for safety reason
        assert ID in _2dCSCG_SpaceAllocator.___space_name___(), \
            " <SpaceInvoker> : space <{}> is not coded yet.".format(ID)
        self._ID_ = ID
        cls_name = _2dCSCG_SpaceAllocator.___space_name___()[ID]
        cls_path = _2dCSCG_SpaceAllocator.___space_path___()[ID]
        self._space_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, inputs, ndim=None, show_info=False):
        if show_info and rAnk == mAster_rank:
            print(f"---[2dCSCG]-[SPACE]-{MyTimer.current_time()}-----")
            print(f"   <space ID>: {self._ID_}")
            print(f"   <space inputs>: {inputs}")

        cOmm.barrier()  # for safety reason
        if ndim is not None: assert ndim == 2
        S = self._space_(inputs, ndim)
        sp = dict()
        sp['type'] = '_2dCSCG_Space'
        sp['ID'] = self._ID_
        sp['inputs'] = inputs
        sp['ndim'] = ndim
        S.___define_parameters___ = sp
        cOmm.barrier()  # for safety reason
        return S





class FormCaller(FrozenOnly):
    """Generate a form."""
    def __init__(self, mesh, space):
        cOmm.barrier()  # for safety reason
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
        cOmm.barrier()  # for safety reason

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
        cOmm.barrier()  # for safety reason
        return FM




class ExactSolutionSelector(FrozenOnly):
    """We select an exact solution object with this class."""
    def __init__(self, mesh):
        cOmm.barrier()  # for safety reason
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, ID, show_info=False, **kwargs):
        if show_info and rAnk == mAster_rank:
            print(f"---[2dCSCG]-[Exact Solution]-{MyTimer.current_time()}-----")
            print(f"   <ES ID>: {ID}")
            print(f"   <ES kwargs>: {kwargs}", flush=True)

        cOmm.barrier()  # for safety reason
        assert ID in _2dCSCG_ExactSolutionAllocator.___exact_solution_name___(), \
            f"Exact solution ID={ID} not found."
        className = _2dCSCG_ExactSolutionAllocator.___exact_solution_name___()[ID]
        classPath = _2dCSCG_ExactSolutionAllocator.___exact_solution_path___()[ID]

        ES =  ExactSolution(self._mesh_)
        status = getattr(import_module(classPath), className)(ES, **kwargs)
        ES.___PRIVATE_set_status___(status)
        esp = dict()
        esp['type'] = '_2dCSCG_ExactSolution'
        esp['ID'] = ID
        esp['mesh_parameters'] = deepcopy(self._mesh_.standard_properties.parameters)
        esp['kwargs'] = kwargs
        ES.___define_parameters___ = esp
        cOmm.barrier()  # for safety reason
        return ES







if __name__ == "__main__":
    # mpiexec -n 4 python objects\CSCG\_2d\master.py

    # mesh = MeshGenerator('cic',)([14,14])
    # mesh = MeshGenerator('rectangle', p_UL=(-1,-1),region_layout=(3,5))([5,5], show_info=True)
    # mesh = MeshGenerator('rectangle_periodic', p_UL=(-1,-1),region_layout=(3,5))([5,5], show_info=True)
    mesh = MeshGenerator('crazy_periodic')([3,3], show_info=True)



    # print(mesh.___PRIVATE_element_division_and_numbering_quality___())
    # mesh.visualize.matplot.element_division()
    #
    mesh.visualize()
    mesh.domain.visualize()
    mesh.domain.regions.visualize()

    # print(rAnk, mesh._element_indices_)

    # space = SpaceInvoker('polynomials')([([-1,0.1,1],), ('Lobatto',1)], show_info=True)
    #
    # # from mifem import save, read
    #
    # mesh.visualize()
    #
    # es = ExactSolutionSelector(mesh)("sL:sincos1", show_info=True)

    # print(mesh.domain.IS.fully_periodic, mesh.domain.boundaries.names)

    # print(MeshGenerator.___domain_input_statistic___())
    # print(MeshGenerator.___domain_input_random_parameters___())