# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class _3dCSCG_DomainPerpendicularSlice(FrozenOnly):
    """"""
    def __init__(self, domain, x, y, z):
        """"""
        regions = domain.regions
        RPS_dict = dict()
        for rn in regions:
            region = regions[rn]

            RSG = region.sub_geometry

            RPS = RSG.make_a_perpendicular_slice_object_on(x=x, y=y, z=z)

            RPS_dict[rn] = RPS

        self._RPS_dict_ = RPS_dict  # will be same in all cores.
        self._freeze_self_()

    @property
    def RPS_dict(self):
        """Region Perpendicular Slice dict.

        A dict whose keys are region names and values are the region Perpendicular Slice. If a
        region has no business with this slice, the value is None.

        """
        return self._RPS_dict_
