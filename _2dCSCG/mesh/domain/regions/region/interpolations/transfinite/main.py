



from _2dCSCG.mesh.domain.regions.region.interpolations.base import InterpolationBase
from _2dCSCG.mesh.domain.regions.region.interpolations.transfinite.mapping import TransfiniteMapping

class Transfinite(InterpolationBase):
    """ The Transfinite interpolation in 2D."""
    def __init__(self, region):
        super().__init__(region)
        assert all([region._edge_geometries_[key].__class__.__name__ != 'Free'
                    for key in region._edge_geometries_]), \
            " <Transfinite> : I do not accept free side geometries. "
        gamma_U = region._edge_geometries_['U'].XY
        gamma_D = region._edge_geometries_['D'].XY
        gamma_L = region._edge_geometries_['L'].XY
        gamma_R = region._edge_geometries_['R'].XY
        dgammaU = region._edge_geometries_['U'].XoYo
        dgammaD = region._edge_geometries_['D'].XoYo
        dgammaL = region._edge_geometries_['L'].XoYo
        dgammaR = region._edge_geometries_['R'].XoYo
        self._TFM_ = TransfiniteMapping(
                (gamma_L, gamma_D, gamma_R, gamma_U),
                (dgammaL, dgammaD, dgammaR, dgammaU))
        self._freeze_self_()

    def mapping(self, r, s):
        """ r, sbe in [0, 1]. """
        r, s = self.___check_rs___(r, s)
        return self._TFM_.mapping(r, s)


    def mapping_X(self, r, s):
        """ r, sbe in [0, 1]. """
        r, s = self.___check_rs___(r, s)
        return self._TFM_.mapping_X(r, s)
    def mapping_Y(self, r, s):
        """ r, sbe in [0, 1]. """
        r, s = self.___check_rs___(r, s)
        return self._TFM_.mapping_Y(r, s)


    def Jacobian_matrix(self, r, s):
        """ r, s be in [0, 1]. """
        r, s = self.___check_rs___(r, s)
        return ((self._TFM_.dx_dr(r, s), self._TFM_.dx_ds(r, s)),
                (self._TFM_.dy_dr(r, s), self._TFM_.dy_ds(r, s)))


    def Jacobian_Xr(self, r, s):
        r, s = self.___check_rs___(r, s)
        return self._TFM_.dx_dr(r, s)
    def Jacobian_Xs(self, r, s):
        r, s = self.___check_rs___(r, s)
        return self._TFM_.dx_ds(r, s)

    def Jacobian_Yr(self, r, s):
        r, s = self.___check_rs___(r, s)
        return self._TFM_.dy_dr(r, s)
    def Jacobian_Ys(self, r, s):
        r, s = self.___check_rs___(r, s)
        return self._TFM_.dy_ds(r, s)