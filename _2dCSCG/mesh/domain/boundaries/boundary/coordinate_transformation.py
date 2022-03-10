


from screws.freeze.main import FrozenOnly


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