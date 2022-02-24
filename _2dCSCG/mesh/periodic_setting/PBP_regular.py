# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from SCREWS.frozen import FrozenOnly



class _2dCSCG_Regular_PBP_RegionSidePair(FrozenOnly):
    def __init__(self, baseMesh, regionEdgeOne, regionEdgeTwo):
        self._baseMesh_ = baseMesh
        rn1, sn1 = regionEdgeOne.split('-')
        assert rn1 in self._baseMesh_.domain.regions.names, \
            " region_name: {} not valid, should be one of <{}>.".format(
                rn1, self._baseMesh_.domain.regions.names)
        assert sn1 in ('U', 'D', 'L', 'R')
        rn2, sn2 = regionEdgeTwo.split('-')
        assert rn2 in self._baseMesh_.domain.regions.names, \
            " region_name: {} not valid, should be one of <{}>.".format(
                rn2, self._baseMesh_.domain.regions.names)
        assert sn2 in ('U', 'D', 'L', 'R')
        assert {'U': 'D', 'D': 'U', 'L': 'R', 'R': 'L'}[sn1] == sn2, \
            " Regular periodic regions side pair must be UD or LR pair!"
        self._regionEdgeOne_ = self._baseMesh_.domain.regions(rn1).edges[sn1]
        self._regionEdgeTwo_ = self._baseMesh_.domain.regions(rn2).edges[sn2]
        self._freeze_self_()

    @property
    def category(self):
        return "regular"

    @property
    def correspondance_of_element_edges(self):
        """
        Returns
        -------
        output: set
            For example:
                {'0U|D2', '3L|R5', ...}
            Note that '|' means regular pairing! Like an internal pairing!

        """
        rn1, sn1 = self._regionEdgeOne_.position.split('-')
        rn2, sn2 = self._regionEdgeTwo_.position.split('-')
        MEGN = self._baseMesh_._element_global_numbering_
        MEGN_R1 = MEGN[rn1]
        MEGN_R2 = MEGN[rn2]
        if sn1 == 'U':
            ElementsSideNumbering1 = MEGN_R1[0, :]
        elif sn1 == 'D':
            ElementsSideNumbering1 = MEGN_R1[-1, :]
        elif sn1 == 'L':
            ElementsSideNumbering1 = MEGN_R1[:, 0]
        elif sn1 == 'R':
            ElementsSideNumbering1 = MEGN_R1[:, -1]
        else:
            raise Exception()
        if sn2 == 'U':
            ElementsSideNumbering2 = MEGN_R2[0, :]
        elif sn2 == 'D':
            ElementsSideNumbering2 = MEGN_R2[-1, :]
        elif sn2 == 'L':
            ElementsSideNumbering2 = MEGN_R2[:, 0]
        elif sn2 == 'R':
            ElementsSideNumbering2 = MEGN_R2[:, -1]
        else:
            raise Exception()
        ESN1 = ElementsSideNumbering1.ravel('F')
        ESN2 = ElementsSideNumbering2.ravel('F')
        CES = list()
        for i, numbering1 in enumerate(ESN1):
            numbering2 = ESN2[i]
            if numbering1 < numbering2:
                CES.append(str(numbering1) + '-' + sn1 + '|' + sn2 + '-' + str(numbering2))
            else:
                CES.append(str(numbering2) + '-' + sn2 + '|' + sn1 + '-' + str(numbering1))
        return CES






class _2dCSCG_Regular_PBP(FrozenOnly):
    def __init__(self, PDS, thePair):
        self._PDS_ = PDS
        self._baseMesh_ = PDS._baseMesh_
        self._boundaryPair_ = thePair
        edgeOne, edgeTwo = thePair.split('=')
        regionEdgesOne =  self._baseMesh_.domain.domain_input.boundary_region_edges[edgeOne]
        regionEdgesTwo =  self._baseMesh_.domain.domain_input.boundary_region_edges[edgeTwo]
        self.___CHECK_REGION_SIDE_PAIRS___(regionEdgesOne, regionEdgesTwo)
        self._freeze_self_()

    def ___CHECK_REGION_SIDE_PAIRS___(self, regionEdgesOne, regionEdgesTwo):
        assert len(regionEdgesOne) == len(regionEdgesTwo), \
            f"different regions sides contained at two sides of regular periodic boundary pair {self._boundaryPair_}."

        self._region_edge_pairs_ = dict()
        for i in range(len(regionEdgesOne)):
            regionEdgeOne = regionEdgesOne[i]
            regionEdgeTwo = regionEdgesTwo[i]
            self._region_edge_pairs_[regionEdgeOne + '=' + regionEdgeTwo] = \
                _2dCSCG_Regular_PBP_RegionSidePair(self._baseMesh_, regionEdgeOne, regionEdgeTwo)

    @property
    def category(self):
        return "regular"

    @property
    def region_edge_pairs(self):
        return self._region_edge_pairs_
