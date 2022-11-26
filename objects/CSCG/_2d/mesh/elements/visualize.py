# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""


from components.freeze.base import FrozenOnly


class _2dCSCG_MeshElements_VIS(FrozenOnly):
    """"""
    def __init__(self, elements):
        self._elements_ = elements
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    def matplot(self, density):
        """To plot the local mesh-elements in the whole computational domain."""
        raise NotImplementedError()
