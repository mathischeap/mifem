# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Tue May 21 17:54:14 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft,
         Delft, the Netherlands

"""
from screws.freeze.main import FrozenOnly
from importlib import import_module


class EdgeGeometryDispatcher(FrozenOnly):
    """ """
    def __init__(self, et):
        """ """
        assert et[0] in self.___coded_edge_geometries___(), \
            " <Region> : EdgeGeometry: {} is not coded. ".format(et[0])
        self._et_ = et
        cls_name = self.___coded_edge_geometries___()[et[0]]
        cls_path = self.___edge_geometries_path___()[et[0]]
        self._cls_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, cc):
        return self._cls_(cc, self._et_)

    @classmethod
    def ___coded_edge_geometries___(cls):
        """ update this whenever we code a new EdgeGeometry. """
        return {'straight': 'Straight',
                'free': 'Free',
                'acw': 'ACW',
                'aacw': 'AACW',
                'customized': 'Customized'}

    @classmethod
    def ___edge_geometries_path___(cls):
        """ """
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        # bast_path = '_2dCSCG.mesh.domain.regions.region.edge_geometries.'
        return {'straight': base_path + 'straight',
                'free': base_path + 'free',
                'acw': base_path + 'acw',
                'aacw': base_path + 'aacw',
                'customized': base_path + 'customized'}