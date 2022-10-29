# -*- coding: utf-8 -*-
"""
This is actually a component of the baseMesh, but we put it here since it is newly
coded for the MPI mesh.

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.periodicSetting.periodic_boundary_pairs.regular.main import _3dCSCG_Regular_PBP



class _3dCSCG_PeriodicDomainSetting(FrozenOnly):
    """ This class is only initialized once in the MasterCore."""
    def __init__(self, baseMesh, givenPairs):
        """ """
        self._baseMesh_ = baseMesh
        self._periodic_boundary_pairs_keys_ = givenPairs
        self.___INITIALIZING_INDIVIDUAL_BOUNDARY_PAIRS___()
        self._periodic_region_side_pairs_ = None
        self._freeze_self_()

    def ___INITIALIZING_INDIVIDUAL_BOUNDARY_PAIRS___(self):
        """ """
        self._periodic_boundary_pairs_ = dict()
        for pair in self._periodic_boundary_pairs_keys_:
            if self._periodic_boundary_pairs_keys_[pair] == 'regular':
                self._periodic_boundary_pairs_[pair] = _3dCSCG_Regular_PBP(self, pair)
            # add more types of periodicBoundaryPair here by adding more elif.
            else:
                raise NotImplementedError(
                    f'Not coded for {self._periodic_boundary_pairs_[pair]} '
                    f'periodic boundary pair.')

    @property
    def periodic_boundary_pairs(self):
        return self._periodic_boundary_pairs_

    @property
    def periodic_region_side_pairs(self):
        if self._periodic_region_side_pairs_ is None:
            self._periodic_region_side_pairs_ = dict()
            for pbp in self.periodic_boundary_pairs:
                self._periodic_region_side_pairs_.update(
                    self.periodic_boundary_pairs[pbp].region_side_pairs)
        return self._periodic_region_side_pairs_