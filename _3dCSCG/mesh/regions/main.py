# -*- coding: utf-8 -*-
"""
A template for Region Class.

@author: Yi Zhang, Created on Mon Sep  3 15:59:16 2018
         Aerodynamics, AE
         TU Delft
"""
from SCREWS.frozen import FrozenOnly
from _3dCSCG.mesh.regions.visualize.main import _3dCSCG_Regions_Visualize


from root.config import *
from typing import Dict, Tuple



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
            the ith side of the regions named `region_name` is on domain boundary.

        """
        return self._domain_._region_sides_on_domain_boundaries_

    @property
    def internal_side_pairs(self):
        """
        Returns
        -------
        self._domain_._region_internal_side_pairs_ : tuple
            A tuple of sets of two elements. Each set represents two internal regions sides which are paired up.

        """
        return self._domain_._region_internal_side_pairs_

    @property
    def map(self):
        """
        Return the regions map.

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
                assert i > 0, f"regions[{now_look_at}] is not connected."

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