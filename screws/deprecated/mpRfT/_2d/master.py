# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/17/2022 9:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from copy import deepcopy
from importlib import import_module
from screws.freeze.main import FrozenOnly
from root.config.main import RANK, MASTER_RANK
from screws.miscellaneous.timer import MyTimer
from objects.CSCG._2d.mesh.domain.inputs.allocator import DomainInputFinder
from objects.CSCG._2d.mesh.domain.main import _2dCSCG_Domain
from objects.CSCG._2d.mesh.main import _2dCSCG_Mesh

from objects.mpRfT._2d.mesh.main import mpRfT2_Mesh

from objects.mpRfT._2d.forms.allocator import mpRfT2_FormsAllocator
from objects.mpRfT._2d.cf.allocator import mpRfT2_FieldsAllocator

from objects.mpRfT._2d.exact_solutions.status.allocator import mpRfT2_ExactSolutionAllocator
from objects.mpRfT._2d.exact_solutions.main import ExactSolution



class MeshGenerator(FrozenOnly):
    def __init__(self, ID, **kwargs):
        """Remember, **kwargs are parameters to customize the domain. The rule is: they can not change the
        topology of the regions!
        """
        di = DomainInputFinder(ID)(**kwargs)
        self._domain_ = _2dCSCG_Domain(di)
        self._ID_ = ID
        self._kwargs_ = kwargs
        self._freeze_self_()

    def __call__(self,
        base_element_layout, dN,
        EDM=None, show_info=False, rfd=None
        ):
        """"""
        if show_info and RANK == MASTER_RANK:
            print(f"---[2dCSCG]-[MESH]-{MyTimer.current_time()}-----")
            print(f"   <domain ID>: {self._ID_}")
            str_dp = str(self._kwargs_)
            if len(str_dp) > 40: str_dp = str_dp[:40] + '...'
            print( "   <domain_parameters>: {}".format(str_dp))
            print(f"   <EDM>: {EDM}", flush=True)

        #------- make a 2d cscg mesh as the base mesh ---------------------------------------
        cscg = _2dCSCG_Mesh(self._domain_, base_element_layout, EDM=EDM)
        mp = dict()
        mp['type'] = '_2dCSCG_Mesh'
        mp['ID'] = self._ID_
        mp['domain_parameters'] = deepcopy(self._kwargs_)
        mp['element_layout'] = base_element_layout
        mp['EDM'] = EDM
        cscg.___define_parameters___ = mp # save the info for restore the cscg mesh

        dp = dict()
        dp['ID'] = self._ID_
        for key in self._kwargs_:
            dp[key] = self._kwargs_[key]
        cscg.domain.___define_parameters___ = dp

        if show_info and RANK == MASTER_RANK:
            str_element_layout = str(base_element_layout)
            if len(str_element_layout) < 40:
                print( "   <cscg element_layout input>: {}".format(str_element_layout))
            else:
                print( "   <cscg element_layout input>: {}...".format(str_element_layout[:40]))
            for rn in cscg.elements.layout:
                print(f"   <cscg element_layout>: {rn} {cscg.elements.layout[rn]}")
            print(f"   <total cscg base elements>: {cscg.elements.GLOBAL_num}", flush=True)

        #------- use the 2d cscg base mesh to make the 2d nCSCG RF2 mesh -----------------------
        mesh = mpRfT2_Mesh(cscg, dN, rfd)

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




class FormCaller(FrozenOnly):
    """"""
    def  __init__(self, mesh):
        self._mesh_ = mesh
        PATH = dict()
        NAME = dict()
        for A in (mpRfT2_FormsAllocator, mpRfT2_FieldsAllocator):
            PATH.update(A.___form_path___())
            NAME.update(A.___form_name___())
        self._PATH_ = PATH
        self._NAME_ = NAME

        self._freeze_self_()

    def __call__(self, ID, *args, **kwargs):
        path = self._PATH_[ID]
        name = self._NAME_[ID]
        CLASS = getattr(import_module(path), name)
        return CLASS(self._mesh_, *args, **kwargs)





class ExactSolutionSelector(FrozenOnly):
    """We select an exact solution object with this class."""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, ID, show_info=False, **kwargs):
        if show_info and RANK == MASTER_RANK:
            print(f"---[mpRfT2]-[Exact Solution]-{MyTimer.current_time()}-----")
            print(f"   <ES ID>: {ID}")
            print(f"   <ES kwargs>: {kwargs}", flush=True)

        assert ID in mpRfT2_ExactSolutionAllocator.___exact_solution_name___(), \
            f"Exact solution ID={ID} not found."
        className = mpRfT2_ExactSolutionAllocator.___exact_solution_name___()[ID]
        classPath = mpRfT2_ExactSolutionAllocator.___exact_solution_path___()[ID]

        ES =  ExactSolution(self._mesh_)
        status = getattr(import_module(classPath), className)(ES, **kwargs)
        ES.___PRIVATE_set_status___(status)
        esp = dict()
        esp['type'] = '_2dCSCG_ExactSolution'
        esp['ID'] = ID
        esp['mesh_parameters'] = deepcopy(self._mesh_.standard_properties.parameters)
        esp['kwargs'] = kwargs
        ES.___define_parameters___ = esp

        return ES



if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/master.py

    mesh = MeshGenerator('rectangle')([3, 3], 2, show_info=False)
    fc = FormCaller(mesh)

    f = fc('0-f-i')

    def p(t, x, y): return x + y + t
    def q(t, x, y): return x + y - t

    s = fc('scalar', p)
    v = fc('vector', (p, q))

    f1i = fc('1-f-i')
    f1o = fc('1-f-o')

    f2i = fc('2-f-i')
    f2o = fc('2-f-o')

    ES = ExactSolutionSelector(mesh)

    es = ES('Poisson:sincos1')

    print(es.status.kinetic_energy(0))

    es.status.kinetic_energy_distribution.visualization()