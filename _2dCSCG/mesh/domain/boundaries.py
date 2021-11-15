# -*- coding: utf-8 -*-
"""

"""
from SCREWS.frozen import FrozenOnly


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



class Boundary(FrozenOnly):
    def __init__(self, boundaries, name):
        """ """
        assert boundaries.ndim == 2, " <Domain> <Boundaries> <Boundary> "
        assert boundaries.__class__.__name__ == '_2dCSCG_Domain_Boundaries'
        assert name in boundaries.names, " <Domain> <Boundaries> <Boundary> "  # meaningless
        self._boundaries_ = boundaries
        self._domain_ = boundaries._domain_
        self._name_ = name
        self._ct_ = BoundaryCoordinateTransformation(self)
        self._freeze_self_()

    @property
    def ndim(self):
        return 2

    @property
    def coordinate_transformation(self):
        """ """
        return self._ct_

    @property
    def position(self):
        return self._boundaries_.region_edges[self.name]

    @property
    def name(self):
        return self._name_



class BoundaryCoordinateTransformation(FrozenOnly):
    def __init__(self, boundary):
        """ """
        assert boundary.ndim == 2, " <Domain> <BoundaryCoordinateTransformation> "
        assert boundary.__class__.__name__ == 'Boundary', " <Domain> <BoundaryCoordinateTransformation> "
        self._boundary_ = boundary
        self._domain_ = boundary._boundaries_._domain_
        self._freeze_self_()

    @property
    def ndim(self):
        return 2