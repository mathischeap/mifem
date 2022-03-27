




from objects.CSCG._2d.mesh.domain.regions.region.edge_geometries.base import EdgeGeometryBase



class Free(EdgeGeometryBase):
    """A free edge geometry is a line we do not really care about its classfication.

    A regions of such edge geometry(s) usually will call a specific interpolator to
    generate the mapping and therefore Jacobian and so on, like the crazy mapping. The
    crazy mapping is a analytical mapping that we only need to know the bounds and c
    which is stored in the `domain_input`.

    While for the transfinite mapping, the edge geometries is essential. We can not set
    it to be free.

    BUT, even for free EdgeGeometry (or SideGeometry in 3D), the corner_coordinates are
    still necessary because they are used to study the topology of the regions.

    Free edge does not mean it is straight! it can be everything (we just donot care).

    """
    def __init__(self, cc, st):
        """ """
        super().__init__(cc, st)
        assert self.edge_type == ('free',), \
            " <EdgeGeometry> <Free> : edge_type={} wrong.".format(self.edge_type)

    def X(self, o):
        raise NotImplementedError()

    def Xo(self, o):
        raise NotImplementedError()

    def Y(self, o):
        raise NotImplementedError()

    def Yo(self, o):
        raise NotImplementedError()
