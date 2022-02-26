# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.frozen import FrozenOnly

from _3dCSCG.mesh.periodic_setting.periodic_boundary_pairs.regular.region_side_pair import _3dCSCG_Regular_PBP_RegionSidePair


class _3dCSCG_Regular_PBP(FrozenOnly):
    def __init__(self, PDS, thePair):
        self._PDS_ = PDS
        self._baseMesh_ = PDS._baseMesh_
        self._boundaryPair_ = thePair
        sideOne, sideTwo = thePair.split('=')
        regionSidesOne =  self._baseMesh_.domain.domain_input.boundary_region_sides[sideOne]
        regionSidesTwo =  self._baseMesh_.domain.domain_input.boundary_region_sides[sideTwo]
        self.___CHECK_REGION_SIDE_PAIRS___(regionSidesOne, regionSidesTwo)
        self._freeze_self_()

    def ___CHECK_REGION_SIDE_PAIRS___(self, regionSidesOne, regionSidesTwo):
        assert len(regionSidesOne) == len(regionSidesTwo), \
            f"different regions sides contained at two sides of regular periodic boundary pair {self._boundaryPair_}."

        self._region_side_pairs_ = dict()
        for i in range(len(regionSidesOne)):
            regionSideOne = regionSidesOne[i]
            regionSideTwo = regionSidesTwo[i]
            self._region_side_pairs_[regionSideOne + '=' + regionSideTwo] = \
                _3dCSCG_Regular_PBP_RegionSidePair(self._baseMesh_, regionSideOne, regionSideTwo)

    @property
    def category(self):
        return "regular"

    @property
    def region_side_pairs(self):
        return self._region_side_pairs_
