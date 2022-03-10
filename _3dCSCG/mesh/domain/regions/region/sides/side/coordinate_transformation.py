






from screws.freeze.main import FrozenOnly


class SideCoordinateTransformation(FrozenOnly):
    def __init__(self, side):
        assert side.ndim == 3, " <Region> <DIII> <SideCoordinateTransformation> "
        assert side.__class__.__name__ == 'Side', " <Region> <DIII> <SideCoordinateTransformation> "
        self._side_ = side
        self._freeze_self_()

    @property
    def ndim(self):
        return 3


    def mapping(self, *evaluation_points):
        """"""
        ep = self._side_.___generate_full_ep___(evaluation_points)
        xyz = self._side_._sides_._region_.interpolation.mapping(*ep)
        return xyz