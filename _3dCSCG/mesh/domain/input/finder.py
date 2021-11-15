# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Fri Dec 14 20:40:12 2018
Aerodynamics, AE
TU Delft
"""
from SCREWS.frozen import FrozenOnly
from SCREWS.exceptions import MeshError
from importlib import import_module



class DomainInputFinder(FrozenOnly):
    """ We use this finder to get a `DomainInput`."""
    def __init__(self, ID):
        try:
            mesh_class = self.___defined_DI___()[ID]
        except KeyError:
            raise MeshError(" <DomainInputFinder> : mesh ID = {} is wrong.".format(ID))
        cls_name = mesh_class
        cls_path = self.___DI_path___()
        self._DomainInput_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()
    
    def __call__(self, *args, **kwargs):
        """"""
        return self._DomainInput_(*args, **kwargs)
    
    @classmethod
    def ___defined_DI___(cls):
        """Here we store all defined meshComponents. Whenever we define a new meshComponents (actually, a new
        domain_input), we add a nickname for it here.
        
        """
        _dict_ = {'crazy': "Crazy",
                  'crazy_periodic': "CrazyPeriodic",
                  'bridge_arch_cracked': "BridgeArchCracked",
                  'psc': "Periodic_Square_Channel",
                  'pwc': "Parallel_Wall_Channel",
                  'LDC': "Lid_Driven_Cavity",
                  }
        return _dict_
        
    @classmethod
    def ___DI_path___(cls):
        """ """
        return '_3dCSCG.mesh.domain.input.domain_inputs'