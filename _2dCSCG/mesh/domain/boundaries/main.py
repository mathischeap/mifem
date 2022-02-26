# -*- coding: utf-8 -*-
"""

"""
from screws.frozen import FrozenOnly
from _2dCSCG.mesh.domain.boundaries.boundary.main import Boundary


class _2dCSCG_Domain_Boundaries(FrozenOnly):
    """ """
    def __init__(self, domain):
        assert domain.ndim == 2, " <Domain> <Boundaries> "
        assert domain.__class__.__name__ == '_2dCSCG_Domain', " <Domain> <Boundaries> "
        self._domain_ = domain
        self._boundaries_ = {}
        for bn in self.names:
            self._boundaries_[bn] = Boundary(self, bn)
        self._freeze_self_()

    def __getitem__(self, bn):
        """ """
        return self._boundaries_[bn]

    @property
    def ndim(self):
        return 2

    @property
    def num(self):
        return self._domain_._num_boundaries_

    @property
    def names(self):
        return self._domain_._boundary_names_

    @property
    def region_edges(self):
        return self._domain_.domain_input.boundary_region_edges



