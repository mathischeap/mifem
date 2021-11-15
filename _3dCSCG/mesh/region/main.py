# -*- coding: utf-8 -*-
"""
A template for Region Class.

@author: Yi Zhang, Created on Mon Sep  3 15:59:16 2018
         Aerodynamics, AE
         TU Delft
"""
from SCREWS.frozen import FrozenOnly
from SCREWS.quadrature import Quadrature
from _3dCSCG.mesh.region.interpolation.searcher import InterpolationSearcher
from _3dCSCG.mesh.region.side_geometry.dispatcher import SideGeometryDispatcher
from _3dCSCG.mesh.region.type_wrt_metric.giver import TypeWr2MetricGiver
from _3dCSCG.mesh.region.sides import Sides
from _3dCSCG.mesh.region.visualize import _3dCSCG_Regions_Visualize

from _3dCSCG.mesh.region.sub_geometry.main import RegionSubGeometry

from root.config import *
from typing import Dict, Tuple


class RegionTopology:
    @classmethod
    def num_sides(cls):
        return 6

    @classmethod
    def num_edges(cls):
        return 12

    @classmethod
    def num_corners(cls):
        return 8

    @classmethod
    def _side_name_to_index_(cls, _side_name_):
        return {'N': 0, 'S': 1, 'W': 2, 'E': 3, 'B': 4, 'F': 5}[_side_name_]

    @classmethod
    def _side_index_to_name_(cls, _side_index_):
        return {0: 'N', 1: 'S', 2: 'W', 3: 'E', 4: 'B', 5: 'F'}[_side_index_]

    @classmethod
    def _side_name_to_index_dict_(cls):
        return {'N': 0, 'S': 1, 'W': 2, 'E': 3, 'B': 4, 'F': 5}

    @classmethod
    def _side_index_to_name_dict_(cls):
        return {0: 'N', 1: 'S', 2: 'W', 3: 'E', 4: 'B', 5: 'F'}

    @classmethod
    def _side_pairing_(cls):
        return {'N': 'S', 'S': 'N', 'W': 'E', 'E': 'W', 'B': 'F', 'F': 'B'}

    @classmethod
    def _side_index_pairing_(cls):
        return {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}

    @classmethod
    def _corner_name_to_index_(cls, _corner_name_):
        return {'NWB': 0, 'SWB': 1, 'NEB': 2, 'SEB': 3, 'NWF': 4, 'SWF': 5, "NEF": 6, 'SEF': 7,
                'NBW': 0, 'SBW': 1, 'NBE': 2, 'SBE': 3, 'NFW': 4, 'SFW': 5, "NFE": 6, 'SFE': 7,
                'WNB': 0, 'WSB': 1, 'ENB': 2, 'ESB': 3, 'WNF': 4, 'WSF': 5, "ENF": 6, 'ESF': 7,
                'WBN': 0, 'WBS': 1, 'EBN': 2, 'EBS': 3, 'WFN': 4, 'WFS': 5, "EFN": 6, 'EFS': 7,
                'BNW': 0, 'BSW': 1, 'BNE': 2, 'BSE': 3, 'FNW': 4, 'FSW': 5, "FNE": 6, 'FSE': 7,
                'BWN': 0, 'BWS': 1, 'BEN': 2, 'BES': 3, 'FWN': 4, 'FWS': 5, "FEN": 6, 'FES': 7}[_corner_name_]

    @classmethod
    def _corner_index_to_name_(cls, _corner_index_):
        return {0: 'NWB', 1: 'SWB', 2: 'NEB', 3: 'SEB', 4: 'NWF', 5: 'SWF', 6: "NEF", 7: 'SEF'}[_corner_index_]

    @classmethod
    def _corner_name_to_index_dict_(cls):
        return {'NWB': 0, 'SWB': 1, 'NEB': 2, 'SEB': 3, 'NWF': 4, 'SWF': 5, "NEF": 6, 'SEF': 7,
                'NBW': 0, 'SBW': 1, 'NBE': 2, 'SBE': 3, 'NFW': 4, 'SFW': 5, "NFE": 6, 'SFE': 7,
                'WNB': 0, 'WSB': 1, 'ENB': 2, 'ESB': 3, 'WNF': 4, 'WSF': 5, "ENF": 6, 'ESF': 7,
                'WBN': 0, 'WBS': 1, 'EBN': 2, 'EBS': 3, 'WFN': 4, 'WFS': 5, "EFN": 6, 'EFS': 7,
                'BNW': 0, 'BSW': 1, 'BNE': 2, 'BSE': 3, 'FNW': 4, 'FSW': 5, "FNE": 6, 'FSE': 7,
                'BWN': 0, 'BWS': 1, 'BEN': 2, 'BES': 3, 'FWN': 4, 'FWS': 5, "FEN": 6, 'FES': 7}

    @classmethod
    def _corner_index_to_name_dict_(cls):
        return {0: 'NWB', 1: 'SWB', 2: 'NEB', 3: 'SEB', 4: 'NWF', 5: 'SWF', 6: "NEF", 7: 'SEF'}

    @classmethod
    def _side_corner_local_numbering_(cls, _side_index_):
        """
        Notice the values are always increasing. The reason is the same with the reason
        for `_side_edge_local_numbering_` or `_side_edge_local_numbering_dict_`.

        """
        return {0: (0, 2, 4, 6),  # the north side has corners locally numbered 0, 2, 4, 6.
                1: (1, 3, 5, 7),  # the south side has corners locally numbered 1, 3, 5, 7.
                2: (0, 1, 4, 5),  # the west side has corners locally numbered 0, 1, 4 5.
                3: (2, 3, 6, 7),  # and so on
                4: (0, 1, 2, 3),
                5: (4, 5, 6, 7)}[_side_index_]

    @classmethod
    def _side_corner_local_numbering_dict_(cls):
        """
        Notice the values are always increasing. The reason is the same with the reason
        for `_side_edge_local_numbering_` or `_side_edge_local_numbering_dict_`.

        """
        return {0: (0, 2, 4, 6),  # the north side has corners locally numbered 0, 2, 4, 6.
                1: (1, 3, 5, 7),  # the south side has corners locally numbered 1, 3, 5, 7.
                2: (0, 1, 4, 5),  # the west side has corners locally numbered 0, 1, 4 5.
                3: (2, 3, 6, 7),  # and so on
                4: (0, 1, 2, 3),
                5: (4, 5, 6, 7)}

    @classmethod
    def _side_edge_local_numbering_(cls, _side_index_):
        """
        Notice the values are always increasing. This is because we are consistently
        using the same convention of numbering local geometries. dx goes first, dy then
        , dz the last.

        """
        return {0: (4, 6, 8, 10),  # the north side has edges locally numbered 4,6,8,10
                1: (5, 7, 9, 11),  # the south side has edges locally numbered 5,7,9,11
                2: (0, 2, 8, 9),  # the west side has edges locally numbered 0,2,8,9
                3: (1, 3, 10, 11),  # and so on
                4: (0, 1, 4, 5),
                5: (2, 3, 6, 7)}[_side_index_]

    @classmethod
    def _side_edge_local_numbering_dict_(cls):
        """
        Notice the values are always increasing. This is because we are consistently
        using the same convention of numbering local geometries. dx goes first, dy then
        , dz the last.

        """
        return {0: (4, 6, 8, 10),  # the north side has edges locally numbered 4,6,8,10
                1: (5, 7, 9, 11),  # the south side has edges locally numbered 5,7,9,11
                2: (0, 2, 8, 9),  # the west side has edges locally numbered 0,2,8,9
                3: (1, 3, 10, 11),  # and so on
                4: (0, 1, 4, 5),
                5: (2, 3, 6, 7)}

    @classmethod
    def _side_axis_distribution_(cls):
        """
        Here 'N':(0,0) means the N side is perpendicular to the 0-axis and at
        the starting side. 'S':(0,-1) means the 'D' side is perpendicular to
        the 0-axis but at the end side, so -1. And 'F':(2,-1) means the F side
        is perpendicular to 2-axis and at the end side.

        And so on.

        """
        return {'N': (0, 0), 'S': (0, -1), 'W': (1, 0), 'E': (1, -1), 'B': (2, 0), 'F': (2, -1)}

    @classmethod
    def _axis_indix_dict_(cls):
        return {'x': 0, 'y': 1, 'z': 2}

    @classmethod
    def _edge_name_to_index_(cls, _name_):
        return \
            {'WB': 0, 'EB': 1, 'WF': 2, 'EF': 3, 'NB': 4, 'SB': 5, 'NF': 6, 'SF': 7, 'NW': 8, 'SW': 9, 'NE': 10,
             'SE': 11,
             'BW': 0, 'BE': 1, 'FW': 2, 'FE': 3, 'BN': 4, 'BS': 5, 'FN': 6, 'FS': 7, 'WN': 8, 'WS': 9, 'EN': 10,
             'ES': 11}[
                _name_]

    @classmethod
    def _edge_index_to_name_(cls, _index_):
        return \
            {0: 'WB', 1: 'EB', 2: 'WF', 3: 'EF', 4: 'NB', 5: 'SB', 6: 'NF', 7: 'SF', 8: 'NW', 9: 'SW', 10: 'NE',
             11: 'SE'}[
                _index_]

    @classmethod
    def _edge_name_to_index_dict_(cls):
        return {'WB': 0, 'EB': 1, 'WF': 2, 'EF': 3, 'NB': 4, 'SB': 5, 'NF': 6, 'SF': 7, 'NW': 8, 'SW': 9, 'NE': 10,
                'SE': 11,
                'BW': 0, 'BE': 1, 'FW': 2, 'FE': 3, 'BN': 4, 'BS': 5, 'FN': 6, 'FS': 7, 'WN': 8, 'WS': 9, 'EN': 10,
                'ES': 11}

    @classmethod
    def _edge_index_to_name_dict_(cls):
        return {0: 'WB', 1: 'EB', 2: 'WF', 3: 'EF', 4: 'NB', 5: 'SB', 6: 'NF', 7: 'SF', 8: 'NW', 9: 'SW', 10: 'NE',
                11: 'SE'}

    @classmethod
    def _corner_local_numbering_(cls):
        """ The local numbering of corners of 3D region. """
        return np.array(range(8)).reshape((2, 2, 2), order='F')

    @classmethod
    def _edge_local_numbering_(cls):
        return (np.arange(4).reshape((1, 2, 2), order='F'),
                np.arange(4).reshape((2, 1, 2), order='F') + 4,
                np.arange(4).reshape((2, 2, 1), order='F') + 8)

    @classmethod
    def _side_local_numbering_(cls):
        return (np.arange(2).reshape((2, 1, 1), order='F'),
                np.arange(2).reshape((1, 2, 1), order='F') + 2,
                np.arange(2).reshape((1, 1, 2), order='F') + 4)


class Region(RegionTopology, FrozenOnly):
    def __init__(self, ndim, name, cc, st, interpolator, domain_input):
        """
        Parameters
        ----------
        ndim : int
            Must to be 3.
        name :
        cc : corner coordinates
        st : side type
        interpolator :
        domain_input :
            The DomainInput, but its classname is not 'DomainInput'. It inherit the
            class `DomainInput` but has personal name.

        """
        assert ndim == 3, " < Region> "
        self._ndim_ = ndim
        self._name_ = name
        self._volume_ = None
        self._domain_input_ = domain_input  # can not remove this line
        self.___PRIVATE_set_corner_coordinates___(cc)
        self.___PRIVATE_set_side_types___(st)
        self.___PRIVATE_parse_side_types___()
        self._type_wrt_metric_ = self.___PRIVATE_PARSE_type_wrt_metric___()
        self._interpolation_ = InterpolationSearcher(interpolator)(self)

        self.___is_orthogonal___ = None
        self._sides_ = Sides(self)
        self._sub_geometry_ = RegionSubGeometry(self)
        self._freeze_self_()

    @property
    def type_wrt_metric(self):
        return self._type_wrt_metric_

    @property
    def sides(self):
        return self._sides_

    @property
    def name(self):
        return self._name_

    @property
    def ndim(self):
        return self._ndim_

    @property
    def interpolation(self):
        return self._interpolation_

    @property
    def volume(self):
        """ Return the volume of this region."""
        if self._volume_ is None:
            p = 20  # we use this many high order numerical quadrature.
            qx, qy, qz, quad_weights = Quadrature([p, p, p], category='Gauss').quad_ndim_ravel
            detJ = self.interpolation.Jacobian(qx, qy, qz)
            self._volume_ = np.sum(detJ * quad_weights) / 8  # divided by 8 because [-1,1]^3 -> [0,1]^3
        return self._volume_

    @property
    def corner_coordinates(self):
        return self._corner_coordinates_

    @property
    def sub_geometry(self):
        return self._sub_geometry_

    def ___PRIVATE_PARSE_type_wrt_metric___(self):
        if self._domain_input_.region_type_wr2_metric is None:
            return TypeWr2MetricGiver('chaotic')(self)
        elif self.name in self._domain_input_.region_type_wr2_metric:
            return TypeWr2MetricGiver(self._domain_input_.region_type_wr2_metric[self.name])(self)
        else:
            return TypeWr2MetricGiver('chaotic')(self)

    def ___PRIVATE_set_side_types___(self, st):
        """
        Parameters
        ----------
        st : dict
            contain the info of non-plane(3D), non-straight(2D) sides.
        """
        _side_types_ = [None for _ in range(self.num_sides())]
        for key in st:
            if key.split('-')[0] == self.name:
                _side_types_[self._side_name_to_index_(key.split('-')[1])] = st[key]
        self._side_types_ = tuple(_side_types_)

    def ___PRIVATE_parse_side_types___(self):
        """
        Here we get the 6(3D) side geometries for further getting the region
        interpolation.

        Attributes
        ----------
        self._side_geometries_ :
            For 3D: N -> S -> W -> E -> B -> F.

        """
        SC_LN = {0: [0, 2, 4, 6],  # the north side has corners locally numbered 0, 2, 4, 6.
                 1: [1, 3, 5, 7],  # the south side has corners locally numbered 1, 3, 5, 7.
                 2: [0, 1, 4, 5],  # the west side has corners locally numbered 0, 1, 4 5.
                 3: [2, 3, 6, 7],  # and so on
                 4: [0, 1, 2, 3],
                 5: [4, 5, 6, 7]}
        _sg_ = dict()
        for i in range(self.num_sides()):
            _st_ = self._side_types_[i]
            if _st_ is None:
                _st_ = ('plane',)
            _cc_ = np.array(self.corner_coordinates)[SC_LN[i]]
            _sg_['NSWEBF'[i]] = SideGeometryDispatcher(_st_)(_cc_)
        self._side_geometries_ = _sg_

    def ___PRIVATE_set_corner_coordinates___(self, cc):
        """
        Parameters
        ----------
        cc : tuple
            For 3D:
            np.shape(cc) must be (8, 3), and cc[i] corresponds to NWB, SWB,
            NEB, SEB, NWF, SWF, NEF, SEF corners, and each item is a tuple
            itself represents (x, y, z) coordinates.

        """
        assert np.shape(cc) == (self.num_corners(), self.ndim), \
            " <Region> : coordinates shape={} wrong.".format(np.shape(cc))
        self._corner_coordinates_ = cc





class Regions(FrozenOnly):
    """A branch of regions that forms a domain. """
    def __init__(self, domain, regions):
        self._domain_ = domain
        self._regions_ = regions
        for key in regions: assert regions[key].__class__.__name__ == 'Region'
        self._orthogonality_ = None
        self._visualize_ = _3dCSCG_Regions_Visualize(self)
        self._freeze_self_()

    def __call__(self, rn=None):
        if rn is None:
            return self._regions_
        else:
            return self._regions_[rn]

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
    def num(self):
        return self._domain_._num_regions_

    @property
    def names(self):
        return self._domain_._region_names_

    @property
    def volumes(self):
        v = {}
        for rn in self.names:
            v[rn] = self(rn).volume
        return v

    @property
    def sides_on_domain_boundaries(self):
        """
        Returns
        -------
        self._domain_._region_sides_on_domain_boundaries_ : dict
            if _region_sides_on_domain_boundaries_[region_name][i] = 1, then we know
            the ith side of the region named `region_name` is on domain boundary.

        """
        return self._domain_._region_sides_on_domain_boundaries_

    @property
    def internal_side_pairs(self):
        """
        Returns
        -------
        self._domain_._region_internal_side_pairs_ : tuple
            A tuple of sets of two elements. Each set represents two internal region sides which are paired up.

        """
        return self._domain_._region_internal_side_pairs_

    @property
    def map(self):
        """
        Return the region map.

        Returns
        -------
        self._domain_._region_map_ : dict

        """
        return self._domain_._region_map_

    @property
    def visualize(self):
        return self._visualize_


    def ___PRIVATE_parse_topology_1___(self):

        if rAnk != mAster_rank: return

        MAP = self.map
        names: Tuple[str] = self.names

        region_coordinates_pool: Dict[str] = dict()
        connection_pool: Dict[str] = dict()

        now_look_at = names[0]
        looked = list()
        we_have_found = list()

        while 1:
            looked.append(now_look_at)

            if len(region_coordinates_pool) == 0:
                region_coordinates_pool[now_look_at] = np.array([0,0,0])

            whats = MAP[now_look_at]

            i = 0
            for w in whats:
                if w[:2] == 'R:':
                    if w not in looked:
                        we_have_found.append(w)
                    i += 1

            if len(names) > 1:
                assert i > 0, f"region[{now_look_at}] is not connected."

            for j, w in enumerate(whats):
                if w[:2] == 'R:':
                    side = 'NSWEBF'[j]
                    if w not in region_coordinates_pool:
                        if side == 'N':
                            region_coordinates_pool[w] = region_coordinates_pool[now_look_at] + np.array([-1, 0, 0])
                        elif side == 'S':
                            region_coordinates_pool[w] = region_coordinates_pool[now_look_at] + np.array([1, 0, 0])
                        elif side == 'W':
                            region_coordinates_pool[w] = region_coordinates_pool[now_look_at] + np.array([0, -1, 0])
                        elif side == 'E':
                            region_coordinates_pool[w] = region_coordinates_pool[now_look_at] + np.array([0, 1, 0])
                        elif side == 'B':
                            region_coordinates_pool[w] = region_coordinates_pool[now_look_at] + np.array([0, 0, -1])
                        elif side == 'F':
                            region_coordinates_pool[w] = region_coordinates_pool[now_look_at] + np.array([0, 0, 1])
                        else:
                            raise Exception()
                        connection_pool[now_look_at + '-' + w] = \
                            np.vstack(
                                [region_coordinates_pool[now_look_at],
                                 region_coordinates_pool[w]]).T

            if len(looked) == len(names):
                break

            while 1:
                to_look = we_have_found.pop()
                if to_look not in looked: # we will for sure find one, otherwise, we should have
                    break

            now_look_at = to_look

        return region_coordinates_pool, connection_pool