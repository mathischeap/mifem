
from screws.frozen import FrozenOnly

from _2dCSCG.mesh.domain.boundaries.boundary.coordinate_transformation import BoundaryCoordinateTransformation


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
