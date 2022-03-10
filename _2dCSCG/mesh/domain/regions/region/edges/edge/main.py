
from screws.freeze.main import FrozenOnly
from _2dCSCG.mesh.domain.regions.region.edges.edge.coordinate_transformation import EdgeCoordinateTransformation


class Edge(FrozenOnly):
    """ """

    def __init__(self, edges, edge_name):
        """ """
        assert edges.ndim == 2, " <Region> <DII> <Edge> "
        assert edges.__class__.__name__ == 'Edges', " <Region> <DII> <Edge> "
        assert edge_name in ('U', 'D', 'L', 'R'), " <Region> <DII> <Edge> "
        self._edges_ = edges
        self._region_ = edges._region_
        self._position_ = edges._region_.name + '-' + edge_name
        self._name_ = edge_name
        self._ct_ = None
        self._freeze_self_()

    @property
    def ndim(self):
        return 2

    @property
    def position(self):
        return self._position_

    @property
    def coordinate_transformation(self):
        """ """
        if self._ct_ is None:
            self._ct_ = EdgeCoordinateTransformation(self)
        return self._ct_


