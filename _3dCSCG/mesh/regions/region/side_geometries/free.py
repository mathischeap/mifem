



from _3dCSCG.mesh.regions.region.side_geometries.base import SideGeometryBase


class Free(SideGeometryBase):
    """
    A free side geometry is a surface we do not really care about its
    classification. A regions of such side geometry(s) usually will call a
    specific interpolator to generate the mapping and therefore Jacobian and
    so on, like the crazy mapping. The crazy mapping is a analytical mapping
    that we only need to know the bounds and c which is stored in the
    `domain_input`.

    While for the transfinite mapping, the side geometries is essential. We can
    not set it to be free.

    """

    def __init__(self, cc, st):
        """ """
        super().__init__(cc, st)
        assert self.side_type == ('free',), \
            " <SideGeometry> <Free> : side_type={} wrong.".format(self.side_type)
