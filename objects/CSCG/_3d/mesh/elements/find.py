# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from screws.freeze.base import FrozenOnly


class _3dCSCG_MeshElements_Find(FrozenOnly):
    """"""

    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()

    def element_at_region_corner(self, region_name, which_corner):
        """

        :param region_name:
        :param which_corner:
        :return:
        """
        mesh = self._elements_._mesh_
        RNS = mesh.domain.regions.names
        assert region_name in RNS, f"region_name = {region_name} is illegal."
        region = mesh.domain.regions[region_name]
        corner_names = region._corner_name_to_index_dict_().keys()
        assert which_corner in corner_names, \
            f"which_corner={which_corner} is illegal!"

        # TODO: to be continued