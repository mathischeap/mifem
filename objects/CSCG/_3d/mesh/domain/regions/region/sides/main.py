# -*- coding: utf-8 -*-




from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.domain.regions.region.sides.side.main import Side


class Sides(FrozenOnly):
    def __init__(self, region):
        assert region.ndim == 3, " <Region> <DIII> <Sides> "
        assert region.__class__.__name__ == 'Region', " <Region> <DIII> <Sides> "
        self._region_ = region
        self._sides_cache_ = {}
        self._freeze_self_()

    @property
    def ndim(self):
        return 3

    def __getitem__(self, side_name):
        """ """
        if side_name in self._sides_cache_:
            pass
        else:
            self._sides_cache_[side_name] = Side(self, side_name)
        return self._sides_cache_[side_name]

    @property
    def types(self):
        return self._region_._side_types_

    @property
    def geometries(self):
        return self._region_._side_geometries_



