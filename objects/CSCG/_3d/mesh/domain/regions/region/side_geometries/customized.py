


from objects.CSCG._3d.mesh.domain.regions.region.side_geometries.base import SideGeometryBase




class Customized(SideGeometryBase):
    """ """

    def __init__(self, cc, st):
        """ """
        super().__init__(cc, st)
        assert self.side_type[0] == 'customized', \
            " <SideGeometryPlane> : side_type[0]={} wrong.".format(self.side_type)
        self._melt_self_()
        self._mapping_ = self.side_type[1]
        self._Jacobian_ = self.side_type[2]
        self._freeze_self_()

    # X __
    def X(self, p, q):
        return self._mapping_[0](p, q)

    def Xp(self, p, q):
        return self._Jacobian_[0][0](p, q)

    def Xq(self, p, q):
        return self._Jacobian_[0][1](p, q)

    # Y __
    def Y(self, p, q):
        return self._mapping_[1](p, q)

    def Yp(self, p, q):
        return self._Jacobian_[1][0](p, q)

    def Yq(self, p, q):
        return self._Jacobian_[1][1](p, q)

    # Z __
    def Z(self, p, q):
        return self._mapping_[2](p, q)

    def Zp(self, p, q):
        return self._Jacobian_[2][0](p, q)

    def Zq(self, p, q):
        return self._Jacobian_[2][1](p, q)
