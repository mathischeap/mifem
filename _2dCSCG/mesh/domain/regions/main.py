# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Tue May 21 14:07:51 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft,
         Delft, Netherlands

"""
import numpy as np
from screws.freeze.main import FrozenOnly
from screws.decorators.accepts import accepts





class Regions(FrozenOnly):
    @accepts('self', '_2dCSCG_Domain', dict)
    def __init__(self, domain, regions):
        self._domain_ = domain
        self._regions_ = regions
        for key in regions: assert regions[key].__class__.__name__ == 'Region'
        self.___generate_global_corner_numbering___()
        self.___generate_region_map___()
        self._freeze_self_()

    def __call__(self, i=None):
        if i is None:
            return self._regions_
        else:
            return self._regions_[i]

    def __getitem__(self, rn):
        return self(rn)

    def __contains__(self, rn):
        return rn in self.names

    def __iter__(self):
        for rn in self.names:
            yield rn

    def __len__(self):
        return len(self.names)

    @property
    def ndim(self):
        return self._domain_.ndim

    @property
    def num(self):
        return self._domain_._num_regions_

    @property
    def names(self):
        return self._domain_._region_names_

    def ___generate_global_corner_numbering___(self):
        """
        Here we parse our regions and try to get the global_corner_numbering
        and something else.

        `_region_corner_coordinates_pool_` and `_global_region_corner_numbering_` are
        not very useful because when two sides at the same location but are on domain
        boundary (like a crack), the global corner numbering will be wrong. So,
        generally, we do not use them.

        Attributes
        ----------
        self._region_corner_coordinates_pool_ : tuple
            A tuple where self._region_corner_pool_[i] indicates the coordinates of the
            regions corner globally numbered as i.
        self._global_region_corner_numbering_ : dict
            A dict whose names indicate regions, and whose value[j].ravel('F')
            indicates the global numbering of [UL, DL, UR, DR][j] corner.

        """
        _rcgn_ = {}
        _corner_pool_ = []
        _current_number_ = 0
        for rn in self.names:  # go through all regions
            _rcgn_[rn] = ()
            for i in range(4):  # go through all corners.
                _corner_ = self(rn).corner_coordinates[i]
                _in_, _index_ = self.___check_duplicated_corner___(_corner_pool_, _corner_)
                if _in_:
                    _rcgn_[rn] += (_index_,)
                else:
                    _rcgn_[rn] += (_current_number_,)
                    _current_number_ += 1
                    _corner_pool_.append(_corner_)
            _rcgn_[rn] = np.array(_rcgn_[rn]).reshape((2, 2), order='F')
        self._corner_coordinates_pool_ = tuple(_corner_pool_)
        self._corner_global_numbering_ = _rcgn_

    @staticmethod
    def ___check_duplicated_corner___(pool, corner):
        """
        Here we check if a corner: corner is already recorded in pool.

        Parameters
        ----------
        pool : list
        corner : tuple

        Returns
        -------
        _in_ : pool
        _index_ : int

        """
        if corner in pool:
            return True, pool.index(corner)
        for i in pool:
            if np.sqrt((i[0] - corner[0]) ** 2 + (i[1] - corner[1]) ** 2) <= 1e-13:
                return True, i
        return False, -1

    def ___generate_region_map___(self):
        """
        Generate the regions map and something else.

        Attributes
        ----------
        self._region_internal_side_pairs_ : tuple
            A tuple of sets of two elements. Each set represents two internal
            regions sides which are paired up.
        self._region_sides_on_domain_boundaries_ : dict
        self._region_map_ : dict

        """
        _rm_ = {}
        # we first find the internal pairing___________________________________________
        for rn in self.names:  # go through all regions
            _rm_[rn] = [[] for _ in range(4)]
            for i in range(4):  # go through all 4 edges of each regions.
                self_corner_indices = self.___found_edge_corner_global_numbering___(rn, i)
                for rnrn in self.names:  # go through all regions except self
                    if rnrn != rn:
                        for ii in range(4):  # go through all 4 edges of the regions.
                            other_corner_indices = self.___found_edge_corner_global_numbering___(rnrn, ii)
                            # noinspection PyTypeChecker
                            if all(self_corner_indices == other_corner_indices):
                                _rm_[rn][i].append(rnrn)
        # We then find the sides on the domain boundaries______________________________
        for bn in self._domain_._boundary_names_:  # we of course go through all domain boundaries
            for db_i in self._domain_.domain_input.boundary_region_edges[bn]:
                _region_name_, _region_edge_ = db_i.split('-')
                _rm_[_region_name_][self(rn)._edge_name_to_index_(_region_edge_)].append(bn)
        # Now we check the regions map and extract more info____________________________
        # We first check each edge only appears at one place and extract the sides on domain boundaries
        _rsodb_ = {}
        for rn in self.names:  # go through all regions
            _rsodb_[rn] = [0 for _ in range(4)]
            for i in range(4):  # go through all 4 edges of each regions.
                try:
                    assert np.shape(_rm_[rn][i]) == (1,)
                except AssertionError:
                    assert np.shape(_rm_[rn][i]) == (2,), \
                        " <Domain>  <2D>: region_map[{}][{}] = {} is wrong, check the domain_input.".format(
                            rn, i, _rm_[rn][i])
                    # noinspection PyUnresolvedReferences
                    rmrni0, rmrni1 = _rm_[rn][i]
                    if rmrni0 in self.names:
                        assert rmrni1 in self._domain_._boundary_names_
                        rmrnib = rmrni1
                        rmrnir = rmrni0
                    elif rmrni0 in self._domain_._boundary_names_:
                        assert rmrni1 in self.names, \
                            " <Domain> : something wrong with boundary_region_edges I think."
                        rmrnib = rmrni0
                        rmrnir = rmrni1
                    else:
                        raise Exception()
                    assert np.shape(_rm_[rmrnir][{0: 1, 1: 0, 2: 3, 3: 2}[i]]) == (2,)
                    assert rn in _rm_[rmrnir][{0: 1, 1: 0, 2: 3, 3: 2}[i]]
                    assert rmrnib in _rm_[rmrnir][{0: 1, 1: 0, 2: 3, 3: 2}[i]]
                    _rm_[rn][i] = [rmrnib, ]
                    _rm_[rmrnir][{0: 1, 1: 0, 2: 3, 3: 2}[i]] = [rmrnib, ]
                _rm_[rn][i] = _rm_[rn][i][0]
                if _rm_[rn][i] in self._domain_._boundary_names_:
                    _rsodb_[rn][i] = 1
                    assert _rm_[rn][i] not in self.names, \
                        " <Regions> <2D> : region_map[{}][{}] = {} is wrong, check the domain_input.".format(
                            rn, i, _rm_[rn][i])
                else:
                    assert _rm_[rn][i] in self.names, \
                        " <Domain3D> : region_map[{}][{}] = {} is wrong, check the domain_input.".format(
                            rn, i, _rm_[rn][i])
            _rsodb_[rn] = tuple(_rsodb_[rn])
            _rm_[rn] = tuple(_rm_[rn])
        # now we check internal sides are correctly paired up and extract the pairing
        # info_________________________________________________________________________
        _region_edge_correct_pairing_ = ({0, 1}, {2, 3})
        _risp_ = ()
        for rn in self.names:  # go through all regions
            for i in range(4):  # go through all 4 edges of each region.
                if not _rsodb_[rn][i]:
                    assert {i, _rm_[_rm_[rn][i]].index(rn)} in _region_edge_correct_pairing_, \
                        " <Domain> <2D> : regions['{}']-side[{}] is paired up with regions['{}']-side[{}]".format(
                            rn, i, _rm_[rn][i], _rm_[_rm_[rn][i]].index(rn))
                    # noinspection PyTypeChecker
                    current_pair = {
                        rn + '-' + self(rn)._edge_index_to_name_(i),
                        _rm_[rn][i] + '-' + self(rn)._edge_index_to_name_(_rm_[_rm_[rn][i]].index(rn))}
                    if current_pair not in _risp_:
                        _risp_ += (current_pair,)

        self._internal_edge_pairs_ = _risp_
        self._edges_on_domain_boundaries_ = _rsodb_
        self._region_map_ = _rm_

    def ___found_edge_corner_global_numbering___(self, region_name, edge_index):
        """
        Parameters
        ----------
        region_name : str
        edge_index : int

        Returns
        -------
        output : ndarray

        """
        return self._corner_global_numbering_[region_name].ravel('F')[
            list(self(region_name)._edge_corner_local_numbering_(edge_index))]

    @property
    def map(self):
        return self._region_map_

    @property
    def internal_edge_pairs(self):
        """
        Returns
        -------
        self._domain_._region_internal_side_pairs_ : tuple
            A tuple of sets of two elements. Each set represents two internal
            regions sides which are paired up.

        """
        return self._internal_edge_pairs_

    @property
    def edges_on_domain_boundaries(self):
        """
        Returns
        -------
        self._edges_on_domain_boundaries_ : dict
            if _edges_on_domain_boundaries_[region_name][i] = 1, then
            we know the ith edge of regions named `region_name` is on domain boundary.

        """
        return self._edges_on_domain_boundaries_