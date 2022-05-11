# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/04 4:06 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from copy import deepcopy
from importlib import import_module
from screws.freeze.main import FrozenOnly
from root.config.main import rAnk, mAster_rank
from screws.miscellaneous.timer import MyTimer
from objects.CSCG._2d.mesh.domain.inputs.allocator import DomainInputFinder
from objects.CSCG._2d.mesh.domain.main import _2dCSCG_Domain
from objects.CSCG._2d.mesh.main import _2dCSCG_Mesh

from objects.nCSCG.rf2._2d.mesh.main import _2nCSCG_RF2_Mesh



class MeshGenerator(FrozenOnly):
    def __init__(self, ID, **kwargs):
        """Remember, **kwargs are parameters to customize the domain.
        The rule is: they can not change the topology of the regions!
        """
        di = DomainInputFinder(ID)(**kwargs)
        self._domain_ = _2dCSCG_Domain(di)
        self._ID_ = ID
        self._kwargs_ = kwargs
        self._freeze_self_()

    def __call__(self, base_element_layout, EDM=None, show_info=False):
        """"""
        if show_info and rAnk == mAster_rank:
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

        if show_info and rAnk == mAster_rank:
            str_element_layout = str(base_element_layout)
            if len(str_element_layout) < 40:
                print( "   <cscg element_layout input>: {}".format(str_element_layout))
            else:
                print( "   <cscg element_layout input>: {}...".format(str_element_layout[:40]))
            for rn in cscg.elements.layout:
                print(f"   <cscg element_layout>: {rn} {cscg.elements.layout[rn]}")
            print(f"   <total cscg base elements>: {cscg.elements.GLOBAL_num}", flush=True)

        #------- use the 2d cscg base mesh to make the 2d nCSCG RF2 mesh ---

        mesh = _2nCSCG_RF2_Mesh(cscg)

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
















if __name__ == "__main__":
    # mpiexec -n 8 python objects/nCSCG/rf2/_2d/master.py
    mesh = MeshGenerator('rectangle')([3,3], show_info=True)
