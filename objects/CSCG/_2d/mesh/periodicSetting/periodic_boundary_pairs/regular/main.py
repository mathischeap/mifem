# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.periodicSetting.periodic_boundary_pairs.regular.region_edge_pair import _2dCSCG_Regular_PBP_RegionEdgePair



class _2dCSCG_Regular_PBP(FrozenOnly):
    def __init__(self, PDS, thePair):
        self._PDS_ = PDS
        self._baseMesh_ = PDS._baseMesh_
        self._boundaryPair_ = thePair
        edgeOne, edgeTwo = thePair.split('=')
        regionEdgesOne =  self._baseMesh_.domain.domain_input.boundary_region_edges[edgeOne]
        regionEdgesTwo =  self._baseMesh_.domain.domain_input.boundary_region_edges[edgeTwo]

        self.___CHECK_REGION_EDGE_PAIRS___(regionEdgesOne, regionEdgesTwo)
        self._freeze_self_()

    def ___CHECK_REGION_EDGE_PAIRS___(self, regionEdgesOne, regionEdgesTwo):
        assert len(regionEdgesOne) == len(regionEdgesTwo), \
            f"different regions sides contained at two sides of regular periodic boundary pair {self._boundaryPair_}."

        self._region_edge_pairs_ = dict()
        for i in range(len(regionEdgesOne)):
            regionEdgeOne = regionEdgesOne[i]
            regionEdgeTwo = regionEdgesTwo[i]
            self._region_edge_pairs_[regionEdgeOne + '=' + regionEdgeTwo] = \
                _2dCSCG_Regular_PBP_RegionEdgePair(self._baseMesh_, regionEdgeOne, regionEdgeTwo)

    @property
    def category(self):
        return "regular"

    @property
    def region_edge_pairs(self):
        return self._region_edge_pairs_
