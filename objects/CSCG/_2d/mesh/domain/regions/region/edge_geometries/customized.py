# -*- coding: utf-8 -*-
from objects.CSCG._2d.mesh.domain.regions.region.edge_geometries.base import EdgeGeometryBase

class Customized(EdgeGeometryBase):
    """ """
    def __init__(self, cc, et):
        """ """
        super().__init__(cc, et)
        assert self.edge_type[0] == 'customized', \
            " <SideGeometryPlane> : edge_type={} wrong.".format(self.edge_type[0])
        self._melt_self_()
        self._mapping_ = self.edge_type[1]
        self._Jacobian_ = self.edge_type[2]
        self._freeze_self_()

    def X(self, o):
        return self._mapping_[0](o)

    def Xo(self, o):
        return self._Jacobian_[0](o)

    def Y(self, o):
        return self._mapping_[1](o)

    def Yo(self, o):
        return self._Jacobian_[1](o)

