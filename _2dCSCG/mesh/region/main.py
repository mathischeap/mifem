# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Tue May 21 14:07:51 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
import numpy as np
from SCREWS.frozen import FrozenOnly
from SCREWS.decorators import accepts
from _2dCSCG.mesh.region.interpolation.searcher import InterpolationSearcher
from _2dCSCG.mesh.region.edge_geometry.dispatcher import EdgeGeometryDispatcher
from _2dCSCG.mesh.region.type_wrt_metric.giver import TypeWr2MetricGiver


from _2dCSCG.mesh.region.edges import Edges



class RegionTopology:
    @classmethod
    def num_edges(cls):
        return 4

    @classmethod
    def num_corners(cls):
        return 4

    @classmethod
    def _edge_name_to_index_dict_(cls):
        return {'U': 0, 'D': 1, 'L': 2, 'R': 3}

    @classmethod
    def _edge_name_to_index_(cls, _edge_name_):
        return {'U': 0, 'D': 1, 'L': 2, 'R': 3}[_edge_name_]

    @classmethod
    def _edge_index_to_name_dict_(cls):
        return {0: 'U', 1: 'D', 2: 'L', 3: 'R'}

    @classmethod
    def _edge_index_to_name_(cls, _edge_index_):
        return {0: 'U', 1: 'D', 2: 'L', 3: 'R'}[_edge_index_]

    @classmethod
    def _edge_pairing_(cls):
        return {'U': 'D', 'D': 'U', 'L': 'R', 'R': 'L'}

    @classmethod
    def _edge_index_pairing_(cls):
        return {0: 1, 1: 0, 2: 3, 3: 2}

    @classmethod
    def _corner_name_to_index_dict_(cls):
        return {'UL': 0, 'DL': 1, 'UR': 2, 'DR': 3}

    @classmethod
    def _corner_name_to_index_(cls, _corner_name_):
        return {'UL': 0, 'DL': 1, 'UR': 2, 'DR': 3,
                'LU': 0, 'LD': 1, 'RU': 2, 'RD': 3}[_corner_name_]

    @classmethod
    def _corner_index_to_name_(cls, _corner_index_):
        return {0: 'UL', 1: 'DL', 2: 'UR', 3: 'DR'}[_corner_index_]

    @classmethod
    def _edge_corner_local_numbering_(cls, _edge_index_):
        """
        Notice the values are always increasing.

        """
        return {0: (0, 2),  # the Upper side has corners locally numbered 0, 2.
                1: (1, 3),  # the Down side has corners locally numbered 1, 3.
                2: (0, 1),  # the Left side has corners locally numbered 0, 1.
                3: (2, 3)}[_edge_index_]

    @classmethod
    def _axis_indix_dict_(cls):
        return {'x': 0, 'y': 1}

    @classmethod
    def _edge_axis_distribution_(cls):
        """
        Here 'U':(0,0) means the U edge is perpendicular to the 0-axis and at the
        starting side. 'D':(0,-1) means the 'D' edge is perpendicular to the 0-axis but
        at the end side, so -1. And 'R':(1,-1) means the 'R' edge is perpendicular to
        1-axis and at the end side. And so on.

        """
        return {'U': (0, 0), 'D': (0, -1), 'L': (1, 0), 'R': (1, -1)}




class Region(RegionTopology, FrozenOnly):
    @accepts('self', int, str)
    def __init__(self, ndim, name, cc, et, interpolator, domain_input):
        """
        Paramters
        ---------
        ndim : int
            Must to be 2.
        name :
        cc : corner coordinates
        et :
        interpolator :
        domain_input :
            The DomainInput, but its classname is not 'DomainInput'. It inherit the
            class `DomainInput` but has personal name.

        """
        assert ndim == 2, " < Region> "
        self._ndim_ = ndim
        self._name_ = name
        self._domain_input_ = domain_input  # can not remove this line
        self.___PRIVATE_set_corner_coordinates___(cc)
        self.___PRIVATE_set_edge_types___(et)
        self.___PRIVATE_parse_edge_types___()
        self._type_wrt_metric_ = self.___PRIVATE_PARSE_type_wrt_metric___()
        self._interpolation_ = InterpolationSearcher(interpolator)(self)
        self._edges_ = Edges(self)
        self._freeze_self_()

    @property
    def edges(self):
        return self._edges_

    @property
    def ndim(self):
        return self._ndim_

    @property
    def name(self):
        return self._name_

    def ___PRIVATE_PARSE_type_wrt_metric___(self):
        if self._domain_input_.region_type_wr2_metric is None:
            return TypeWr2MetricGiver('chaotic')(self)
        elif self.name in self._domain_input_.region_type_wr2_metric:
            return TypeWr2MetricGiver(self._domain_input_.region_type_wr2_metric[self.name])(self)
        else:
            return TypeWr2MetricGiver('chaotic')(self)


    @property
    def type_wrt_metric(self):
        return self._type_wrt_metric_

    def ___PRIVATE_set_corner_coordinates___(self, cc):
        """
        Parameters
        ----------
        cc : tuple
            For 2D:
            np.shape(cc) must be (4, 2), and cc[i] corresponds to UL, DL, UR,
            DR corners, and each item is a tuple itself represents (x, y)
            coordinates.

        """
        assert np.shape(cc) == (4, 2), \
            " <Region> : coordinates shape={} wrong.".format(np.shape(cc))
        self._corner_coordinates_ = cc

    @property
    def corner_coordinates(self):
        return self._corner_coordinates_

    @accepts('self', dict)
    def ___PRIVATE_set_edge_types___(self, et):
        """
        Parameters
        ----------
        et : dict
            contain the info of non-straight edges.

        """
        _edge_types_ = [None for _ in range(self.num_edges())]
        for key in et:
            if key.split('-')[0] == self.name:
                _edge_types_[self._edge_name_to_index_(key.split('-')[1])] = et[key]
        self._edge_types_ = tuple(_edge_types_)

    def ___PRIVATE_parse_edge_types___(self):
        """
        Here we get the 4 edge geometries for futher getting the regions interpolation.

        Attributes
        ----------
        self._edge_geometries_ :
            For 2D: U -> D -> L -> R

        """
        _eg_ = {}
        for i in range(self.num_edges()):
            _et_ = self._edge_types_[i]
            if _et_ is None: # when we do not mention, it is straight, not free.
                _et_ = ('straight',)
            _cc_ = np.array(self.corner_coordinates)[list(self._edge_corner_local_numbering_(i))]
            _eg_[self._edge_index_to_name_(i)] = EdgeGeometryDispatcher(_et_)(_cc_)
        self._edge_geometries_ = _eg_

    @property
    def interpolation(self):
        return self._interpolation_

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
            A dict whose names indicate regions, and whose valune[j].ravel('F')
            indicates the golobal numbering of [UL, DL, UR, DR][j] corner.

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
            for i in range(4):  # go through all 4 edges of each regions.
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