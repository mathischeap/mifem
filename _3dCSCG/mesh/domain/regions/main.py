# -*- coding: utf-8 -*-
"""
A template for Region Class.

@author: Yi Zhang, Created on Mon Sep  3 15:59:16 2018
         Aerodynamics, AE
         TU Delft
"""
from screws.freeze.main import FrozenOnly
from _3dCSCG.mesh.domain.regions.visualize.main import _3dCSCG_Regions_Visualize
from _3dCSCG.mesh.domain.regions.topology import _3dCSCG_Regions_Topology





class Regions(FrozenOnly):
    """A branch of regions that forms a domain. """
    def __init__(self, domain, regions):
        self._domain_ = domain
        self._regions_ = regions
        self._topology_ = None
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

    @property
    def topology(self):
        if self._topology_ is None:
            self._topology_ = _3dCSCG_Regions_Topology(self)
        return self._topology_
