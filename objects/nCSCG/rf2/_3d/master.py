# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/04 3:51 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from copy import deepcopy
from importlib import import_module
from screws.freeze.base import FrozenOnly
from screws.miscellaneous.timer import MyTimer
from root.config.main import rAnk, mAster_rank
from objects.CSCG._3d.mesh.domain.inputs.allocator import DomainInputAllocator
from objects.CSCG._3d.mesh.domain.main import _3dCSCG_Domain
from objects.CSCG._3d.mesh.main import _3dCSCG_Mesh

from objects.nCSCG.rf2._3d.mesh.main import _3nCSCG_RF2_Mesh





class MeshGenerator(FrozenOnly):
    def __init__(self, ID, **kwargs):
        """Initialize a mesh generator.

        Remember, **kwargs are parameters to customize the domain.
        The rule is: they can not change the topology of the regions!
        """
        di = DomainInputAllocator(ID)(**kwargs)
        self._domain_ = _3dCSCG_Domain(di)
        self._ID_ = ID
        self._kwargs_ = kwargs
        self._freeze_self_()

    def __call__(self, base_element_layout, EDM=None, show_info=False):
        """Call to make the mesh instance.

        :param base_element_layout:
        :param EDM: ``Element Distribution Method``. When
            EDM = 'debug', we use a naive method
            EDM = None, we find a proper one.
            EDM = a string, we use the specific one.
        :param show_info:
        :return:
        """
        if show_info and rAnk == mAster_rank:
            print(f"---[3dCSCG]-[MESH]-{MyTimer.current_time()}-----")
            print(f"   <domain ID>: {self._ID_}")
            str_dp = str(self._kwargs_)
            if len(str_dp) > 40: str_dp = str_dp[:40] + '...'
            print( "   <domain_parameters>: {}".format(str_dp))
            print(f"   <EDM>: {EDM}", flush=True)

        #-------------------- make a cscg base mesh ------------------------------------------------
        cscg = _3dCSCG_Mesh(self._domain_, base_element_layout, EDM=EDM)

        mp = dict()
        mp['type'] = '_3dCSCG_Mesh'
        mp['ID'] = self._ID_
        mp['domain_parameters'] = deepcopy(self._kwargs_)
        mp['element_layout'] = base_element_layout
        mp['EDM'] = EDM
        cscg.___define_parameters___ = mp # use these parameters to restore a cscg mesh

        dp = dict()
        dp['ID'] = self._ID_
        for key in self._kwargs_:
            dp[key] = self._kwargs_[key]
        cscg.domain.___define_parameters___ = dp

        if show_info and rAnk == mAster_rank:
            str_element_layout = str(base_element_layout)
            if len(str_element_layout) < 40:
                print( f"   <cscg element_layout input>: {str_element_layout}")
            else:
                print( f"   <cscg element_layout input>: {str_element_layout[:40]}...")
            for rn in cscg.elements.layout:
                print(f"   <cscg element_layout>: {rn} {cscg.elements.layout[rn]}")
            print(f"   <total base cscg elements>: {cscg.elements.GLOBAL_num}", flush=True)

        #----------- use the cscg base mesh to make a _3nCSCG_RF2_Mesh -----------------------------
        mesh = _3nCSCG_RF2_Mesh(cscg)

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






if __name__ == "__main__":
    # mpiexec -n 8 python objects/nCSCG/rf2/_3d/master.py
    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None, show_info=True)
