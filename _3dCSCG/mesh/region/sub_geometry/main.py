"""
Sub-geometries of a region.

For example, we can pick a slice from a region. Then we can do reconstruction of a form on this slice to view some
particular structure.

"""


from _3dCSCG.mesh.region.sub_geometry.perpendicular_slice import RegionPerpendicularSlice


from SCREWS.frozen import FrozenOnly


class RegionSubGeometry(FrozenOnly):
    """"""

    def __init__(self, region):
        """"""
        self._region_ = region

        self._freeze_self_()





    def GENERATE_perpendicular_slice_object(self, *args, **kwargs):
        """"""
        return RegionPerpendicularSlice(self._region_, *args, **kwargs)

