# -*- coding: utf-8 -*-


from screws.freeze.main import FrozenOnly


class _3dCSCG_Regular_PBP_RegionSidePair(FrozenOnly):
    def __init__(self, baseMesh, regionSideOne, regionSideTwo):
        self._baseMesh_ = baseMesh
        rn1, sn1 = regionSideOne.split('-')
        assert rn1 in self._baseMesh_.domain.regions.names, \
            " region_name: {} not valid, should be one of <{}>.".format(
                rn1, self._baseMesh_.domain.regions.names)
        assert sn1 in ('N', 'S', 'W', 'E', 'B', 'F')
        rn2, sn2 = regionSideTwo.split('-')
        assert rn2 in self._baseMesh_.domain.regions.names, \
            " region_name: {} not valid, should be one of <{}>.".format(
                rn2, self._baseMesh_.domain.regions.names)
        assert sn2 in ('N', 'S', 'W', 'E', 'B', 'F')
        assert {'N': 'S', 'S': 'N', 'W': 'E', 'E': 'W', 'B': 'F', 'F': 'B'}[sn1] == sn2, \
            " Regular periodic regions side pair must be NS, WE or BF pair!"
        self._regionSideOne_ = self._baseMesh_.domain.regions(rn1).sides[sn1]
        self._regionSideTwo_ = self._baseMesh_.domain.regions(rn2).sides[sn2]
        self._freeze_self_()

    @property
    def category(self):
        return "regular"

    @property
    def correspondence_of_element_sides(self):
        """
        Returns
        -------
        output: set
            For example:
                {'0N|S2', '3N|S5', '6N|S8', '9N|S11', '12N|S14', '15N|S17', '18N|S20
                 ', '21N|S23', '24N|S26'}
            Note that '|' means regular pairing! Like an internal pairing!

        """
        rn1, sn1 = self._regionSideOne_.position.split('-')
        rn2, sn2 = self._regionSideTwo_.position.split('-')
        MEGN = self._baseMesh_._element_global_numbering_
        MEGN_R1 = MEGN[rn1]
        MEGN_R2 = MEGN[rn2]
        if sn1 == 'N':
            ElementsSideNumbering1 = MEGN_R1[0, :, :]
        elif sn1 == 'S':
            ElementsSideNumbering1 = MEGN_R1[-1, :, :]
        elif sn1 == 'W':
            ElementsSideNumbering1 = MEGN_R1[:, 0, :]
        elif sn1 == 'E':
            ElementsSideNumbering1 = MEGN_R1[:, -1, :]
        elif sn1 == 'B':
            ElementsSideNumbering1 = MEGN_R1[:, :, 0]
        elif sn1 == 'F':
            ElementsSideNumbering1 = MEGN_R1[:, :, -1]
        else:
            raise Exception()
        if sn2 == 'N':
            ElementsSideNumbering2 = MEGN_R2[0, :, :]
        elif sn2 == 'S':
            ElementsSideNumbering2 = MEGN_R2[-1, :, :]
        elif sn2 == 'W':
            ElementsSideNumbering2 = MEGN_R2[:, 0, :]
        elif sn2 == 'E':
            ElementsSideNumbering2 = MEGN_R2[:, -1, :]
        elif sn2 == 'B':
            ElementsSideNumbering2 = MEGN_R2[:, :, 0]
        elif sn2 == 'F':
            ElementsSideNumbering2 = MEGN_R2[:, :, -1]
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
