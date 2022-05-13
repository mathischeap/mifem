# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly
from screws.decorators.accepts import memoize1
from objects.CSCG._2d.mesh.domain.regions.region.edges.edge.main import Edge

class Edges(FrozenOnly):
    """ """

    def __init__(self, region):
        """ """
        assert region.ndim == 2, " <Region> <DII> <Edges> "
        assert region.__class__.__name__ == 'Region', " <Region> <DII> <Edges> "
        self._region_ = region
        self._freeze_self_()

    @property
    def ndim(self):
        return 2

    @memoize1  # this decorator is great!
    def __call__(self, edge_name):
        """ """
        return Edge(self, edge_name)

    def __getitem__(self, item):
        """ """
        return self(item)

    @property
    def types(self):
        return self._region_._edge_types_

    @property
    def geometries(self):
        return self._region_._edge_geometries_



