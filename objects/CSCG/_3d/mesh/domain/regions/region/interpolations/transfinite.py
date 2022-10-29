




from root.config.main import *

import random
from objects.CSCG._3d.mesh.domain.regions.region.interpolations.base import InterpolationBase

from screws.exceptions import ThreeDimensionalTransfiniteInterpolationError

from screws.numerical._2dSpace.Jacobian_21 import NumericalPartialDerivative_xy



class Transfinite(InterpolationBase):
    """ The Transfinite interpolation in 3D. """

    def __init__(self, region):
        """
        To initialize a transfinite interpolation, we take a regions as input.

        Parameters
        ----------
        region : Region

        """
        super().__init__(region)
        assert all([region._side_geometries_[key].__class__.__name__ != 'Free'
                    for key in region._side_geometries_]), \
            " <Transfinite> : I do not accept `Free` side geometries. "
        self.___check_line_consistency___()
        self.___check_Jacobian___()
        self._freeze_self_()

    def ___check_line_consistency___(self):
        """
        We know that trace to a line from two sides should result in the same
        line function. So we do this check to see if the provided
        _side_geometries_ are correct.

        Notice that we do not have to check the corner consistency anymore if
        this check is all right.

        """
        I = random.randint(2, 4)
        rst = np.linspace(0, 1, I)
        _NW_, _WN_ = self.___NW___(rst), self.___WN___(rst)
        _NE_, _EN_ = self.___NE___(rst), self.___EN___(rst)
        _NB_, _BN_ = self.___NB___(rst), self.___BN___(rst)
        _NF_, _FN_ = self.___NF___(rst), self.___FN___(rst)
        _SW_, _WS_ = self.___SW___(rst), self.___WS___(rst)
        _SE_, _ES_ = self.___SE___(rst), self.___ES___(rst)
        _SB_, _BS_ = self.___SB___(rst), self.___BS___(rst)
        _SF_, _FS_ = self.___SF___(rst), self.___FS___(rst)
        _WB_, _BW_ = self.___WB___(rst), self.___BW___(rst)
        _WF_, _FW_ = self.___WF___(rst), self.___FW___(rst)
        _EB_, _BE_ = self.___EB___(rst), self.___BE___(rst)
        _EF_, _FE_ = self.___EF___(rst), self.___FE___(rst)

        try:
            for i in range(3):
                np.testing.assert_array_almost_equal(_NW_[i], _WN_[i])
                np.testing.assert_array_almost_equal(_NE_[i], _EN_[i])
                np.testing.assert_array_almost_equal(_NB_[i], _BN_[i])
                np.testing.assert_array_almost_equal(_NF_[i], _FN_[i])
                np.testing.assert_array_almost_equal(_SW_[i], _WS_[i])
                np.testing.assert_array_almost_equal(_SE_[i], _ES_[i])
                np.testing.assert_array_almost_equal(_SB_[i], _BS_[i])
                np.testing.assert_array_almost_equal(_SF_[i], _FS_[i])
                np.testing.assert_array_almost_equal(_WB_[i], _BW_[i])
                np.testing.assert_array_almost_equal(_WF_[i], _FW_[i])
                np.testing.assert_array_almost_equal(_EB_[i], _BE_[i])
                np.testing.assert_array_almost_equal(_EF_[i], _FE_[i])
            raise_Error = False
        except AssertionError:
            raise_Error = True

        raise_Error = COMM.allreduce(raise_Error, op=MPI.LOR)

        if raise_Error:
            raise ThreeDimensionalTransfiniteInterpolationError(
                "Something is wrong with 3d TransfiniteInterpolation")

    def ___check_Jacobian___(self):
        """ """
        p = np.random.rand(7, 5)
        q = np.random.rand(7, 5)
        for side in ('N', 'S', 'W', 'E', 'B', 'F'):
            X = self._region_._side_geometries_[side].X
            Y = self._region_._side_geometries_[side].Y
            Z = self._region_._side_geometries_[side].Z
            Xp = self._region_._side_geometries_[side].Xp
            Xq = self._region_._side_geometries_[side].Xq
            Yp = self._region_._side_geometries_[side].Yp
            Yq = self._region_._side_geometries_[side].Yq
            Zp = self._region_._side_geometries_[side].Zp
            Zq = self._region_._side_geometries_[side].Zq
            XNJ = NumericalPartialDerivative_xy(X, p, q)
            assert all(XNJ.check_total(Xp, Xq)), " <Transfinite> : X, Xp, Xq wrong"
            YNJ = NumericalPartialDerivative_xy(Y, p, q)
            assert all(YNJ.check_total(Yp, Yq)), " <Transfinite> : Y, Yp, Yq wrong"
            ZNJ = NumericalPartialDerivative_xy(Z, p, q)
            assert all(ZNJ.check_total(Zp, Zq)), " <Transfinite> : Z, Zp, Zq wrong"



    def mapping(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        N = self.___NnN___(s, t)
        S = self.___SsS___(s, t)
        W = self.___WwW___(r, t)
        E = self.___EeE___(r, t)
        B = self.___BbB___(r, s)
        F = self.___FfF___(r, s)

        NW = self.___NW___(t)
        NE = self.___NE___(t)
        SW = self.___SW___(t)
        SE = self.___SE___(t)
        NB = self.___NB___(s)
        NF = self.___NF___(s)
        SB = self.___SB___(s)
        SF = self.___SF___(s)
        WB = self.___WB___(r)
        WF = self.___WF___(r)
        EB = self.___EB___(r)
        EF = self.___EF___(r)

        NWB = self.___N_W_B___
        NWF = self.___N_W_F___
        NEB = self.___N_E_B___
        NEF = self.___N_E_F___
        SWB = self.___S_W_B___
        SWF = self.___S_W_F___
        SEB = self.___S_E_B___
        SEF = self.___S_E_F___

        x1 = (1 - r) * N + r * S
        x2 = (1 - s) * W + s * E
        x3 = (1 - t) * B + t * F

        x12 = (1 - r) * (1 - s) * NW + (1 - r) * s * NE + r * (1 - s) * SW + r * s * SE
        x13 = (1 - r) * (1 - t) * NB + (1 - r) * t * NF + r * (1 - t) * SB + r * t * SF
        x23 = (1 - s) * (1 - t) * WB + (1 - s) * t * WF + s * (1 - t) * EB + s * t * EF

        x123 = [None, None, None]
        for i in range(3):
            x123[i] = + (1 - r) * (1 - s) * (1 - t) * NWB[i] \
                      + (1 - r) * (1 - s) * t * NWF[i] \
                      + (1 - r) * s * (1 - t) * NEB[i] \
                      + (1 - r) * s * t * NEF[i] \
                      + r * (1 - s) * (1 - t) * SWB[i] \
                      + r * (1 - s) * t * SWF[i] \
                      + r * s * (1 - t) * SEB[i] \
                      + r * s * t * SEF[i]

        X = x1 + x2 + x3 - x12 - x13 - x23 + x123
        return X

    def mapping_X(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        N = self.___NnN___(s, t)[0]
        S = self.___SsS___(s, t)[0]
        W = self.___WwW___(r, t)[0]
        E = self.___EeE___(r, t)[0]
        B = self.___BbB___(r, s)[0]
        F = self.___FfF___(r, s)[0]

        NW = self.___NW___(t)[0]
        NE = self.___NE___(t)[0]
        SW = self.___SW___(t)[0]
        SE = self.___SE___(t)[0]
        NB = self.___NB___(s)[0]
        NF = self.___NF___(s)[0]
        SB = self.___SB___(s)[0]
        SF = self.___SF___(s)[0]
        WB = self.___WB___(r)[0]
        WF = self.___WF___(r)[0]
        EB = self.___EB___(r)[0]
        EF = self.___EF___(r)[0]

        NWB = self.___N_W_B___[0]
        NWF = self.___N_W_F___[0]
        NEB = self.___N_E_B___[0]
        NEF = self.___N_E_F___[0]
        SWB = self.___S_W_B___[0]
        SWF = self.___S_W_F___[0]
        SEB = self.___S_E_B___[0]
        SEF = self.___S_E_F___[0]

        x1 = (1 - r) * N + r * S
        x2 = (1 - s) * W + s * E
        x3 = (1 - t) * B + t * F

        x12 = (1 - r) * (1 - s) * NW + (1 - r) * s * NE + r * (1 - s) * SW + r * s * SE
        x13 = (1 - r) * (1 - t) * NB + (1 - r) * t * NF + r * (1 - t) * SB + r * t * SF
        x23 = (1 - s) * (1 - t) * WB + (1 - s) * t * WF + s * (1 - t) * EB + s * t * EF

        x123 = + (1 - r) * (1 - s) * (1 - t) * NWB \
               + (1 - r) * (1 - s) * t * NWF \
               + (1 - r) * s * (1 - t) * NEB \
               + (1 - r) * s * t * NEF \
               + r * (1 - s) * (1 - t) * SWB \
               + r * (1 - s) * t * SWF \
               + r * s * (1 - t) * SEB \
               + r * s * t * SEF

        X = x1 + x2 + x3 - x12 - x13 - x23 + x123
        return X

    def mapping_Y(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        N = self.___NnN___(s, t)[1]
        S = self.___SsS___(s, t)[1]
        W = self.___WwW___(r, t)[1]
        E = self.___EeE___(r, t)[1]
        B = self.___BbB___(r, s)[1]
        F = self.___FfF___(r, s)[1]

        NW = self.___NW___(t)[1]
        NE = self.___NE___(t)[1]
        SW = self.___SW___(t)[1]
        SE = self.___SE___(t)[1]
        NB = self.___NB___(s)[1]
        NF = self.___NF___(s)[1]
        SB = self.___SB___(s)[1]
        SF = self.___SF___(s)[1]
        WB = self.___WB___(r)[1]
        WF = self.___WF___(r)[1]
        EB = self.___EB___(r)[1]
        EF = self.___EF___(r)[1]

        NWB = self.___N_W_B___[1]
        NWF = self.___N_W_F___[1]
        NEB = self.___N_E_B___[1]
        NEF = self.___N_E_F___[1]
        SWB = self.___S_W_B___[1]
        SWF = self.___S_W_F___[1]
        SEB = self.___S_E_B___[1]
        SEF = self.___S_E_F___[1]

        x1 = (1 - r) * N + r * S
        x2 = (1 - s) * W + s * E
        x3 = (1 - t) * B + t * F

        x12 = (1 - r) * (1 - s) * NW + (1 - r) * s * NE + r * (1 - s) * SW + r * s * SE
        x13 = (1 - r) * (1 - t) * NB + (1 - r) * t * NF + r * (1 - t) * SB + r * t * SF
        x23 = (1 - s) * (1 - t) * WB + (1 - s) * t * WF + s * (1 - t) * EB + s * t * EF

        x123 = + (1 - r) * (1 - s) * (1 - t) * NWB \
               + (1 - r) * (1 - s) * t * NWF \
               + (1 - r) * s * (1 - t) * NEB \
               + (1 - r) * s * t * NEF \
               + r * (1 - s) * (1 - t) * SWB \
               + r * (1 - s) * t * SWF \
               + r * s * (1 - t) * SEB \
               + r * s * t * SEF

        Y = x1 + x2 + x3 - x12 - x13 - x23 + x123
        return Y

    def mapping_Z(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        N = self.___NnN___(s, t)[2]
        S = self.___SsS___(s, t)[2]
        W = self.___WwW___(r, t)[2]
        E = self.___EeE___(r, t)[2]
        B = self.___BbB___(r, s)[2]
        F = self.___FfF___(r, s)[2]

        NW = self.___NW___(t)[2]
        NE = self.___NE___(t)[2]
        SW = self.___SW___(t)[2]
        SE = self.___SE___(t)[2]
        NB = self.___NB___(s)[2]
        NF = self.___NF___(s)[2]
        SB = self.___SB___(s)[2]
        SF = self.___SF___(s)[2]
        WB = self.___WB___(r)[2]
        WF = self.___WF___(r)[2]
        EB = self.___EB___(r)[2]
        EF = self.___EF___(r)[2]

        NWB = self.___N_W_B___[2]
        NWF = self.___N_W_F___[2]
        NEB = self.___N_E_B___[2]
        NEF = self.___N_E_F___[2]
        SWB = self.___S_W_B___[2]
        SWF = self.___S_W_F___[2]
        SEB = self.___S_E_B___[2]
        SEF = self.___S_E_F___[2]

        x1 = (1 - r) * N + r * S
        x2 = (1 - s) * W + s * E
        x3 = (1 - t) * B + t * F

        x12 = (1 - r) * (1 - s) * NW + (1 - r) * s * NE + r * (1 - s) * SW + r * s * SE
        x13 = (1 - r) * (1 - t) * NB + (1 - r) * t * NF + r * (1 - t) * SB + r * t * SF
        x23 = (1 - s) * (1 - t) * WB + (1 - s) * t * WF + s * (1 - t) * EB + s * t * EF

        x123 = + (1 - r) * (1 - s) * (1 - t) * NWB \
               + (1 - r) * (1 - s) * t * NWF \
               + (1 - r) * s * (1 - t) * NEB \
               + (1 - r) * s * t * NEF \
               + r * (1 - s) * (1 - t) * SWB \
               + r * (1 - s) * t * SWF \
               + r * s * (1 - t) * SEB \
               + r * s * t * SEF

        Z = x1 + x2 + x3 - x12 - x13 - x23 + x123
        return Z

    def Jacobian_matrix(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        N = self.___NnN___(s, t)
        S = self.___SsS___(s, t)
        Ns = self.___Ns_Ns___(s, t)
        Ss = self.___Ss_Ss___(s, t)
        Nt = self.___Nt_Nt___(s, t)
        St = self.___St_St___(s, t)

        W = self.___WwW___(r, t)
        E = self.___EeE___(r, t)
        Wr = self.___Wr_Wr___(r, t)
        Er = self.___Er_Er___(r, t)
        Wt = self.___Wt_Wt___(r, t)
        Et = self.___Et_Et___(r, t)

        B = self.___BbB___(r, s)
        F = self.___FfF___(r, s)
        Br = self.___Br_Br___(r, s)
        Fr = self.___Fr_Fr___(r, s)
        Bs = self.___Bs_Bs___(r, s)
        Fs = self.___Fs_Fs___(r, s)

        NW = self.___NW___(t)
        NE = self.___NE___(t)
        SW = self.___SW___(t)
        SE = self.___SE___(t)
        NWt = self.___NW_t___(t)
        NEt = self.___NE_t___(t)
        SWt = self.___SW_t___(t)
        SEt = self.___SE_t___(t)

        NB = self.___NB___(s)
        NF = self.___NF___(s)
        SB = self.___SB___(s)
        SF = self.___SF___(s)
        NBs = self.___NB_s___(s)
        NFs = self.___NF_s___(s)
        SBs = self.___SB_s___(s)
        SFs = self.___SF_s___(s)

        WB = self.___WB___(r)
        WF = self.___WF___(r)
        EB = self.___EB___(r)
        EF = self.___EF___(r)
        WBr = self.___WB_r___(r)
        WFr = self.___WF_r___(r)
        EBr = self.___EB_r___(r)
        EFr = self.___EF_r___(r)

        NWB = self.___N_W_B___
        NWF = self.___N_W_F___
        NEB = self.___N_E_B___
        NEF = self.___N_E_F___
        SWB = self.___S_W_B___
        SWF = self.___S_W_F___
        SEB = self.___S_E_B___
        SEF = self.___S_E_F___

        x1 = -_1_ * N + _1_ * S
        x2 = (1 - s) * Wr + s * Er
        x3 = (1 - t) * Br + t * Fr
        x12 = -s * NE - (1 - s) * NW + (1 - s) * SW + s * SE
        x13 = -t * NF - (1 - t) * NB + (1 - t) * SB + t * SF
        x23 = (1 - s) * (1 - t) * WBr + (1 - s) * t * WFr + s * (1 - t) * EBr + s * t * EFr
        x123 = [None, None, None]
        for i in range(3):
            x123[i] = - (1 - s) * (1 - t) * NWB[i] \
                      - (1 - s) * t * NWF[i] \
                      - s * (1 - t) * NEB[i] \
                      - s * t * NEF[i] \
                      + (1 - s) * (1 - t) * SWB[i] \
                      + (1 - s) * t * SWF[i] \
                      + s * (1 - t) * SEB[i] \
                      + s * t * SEF[i]
        xr, yr, zr = x1 + x2 + x3 - x12 - x13 - x23 + x123

        x1 = (1 - r) * Ns + r * Ss
        x2 = -_1_ * W + _1_ * E
        x3 = (1 - t) * Bs + t * Fs
        x12 = -(1 - r) * NW + (1 - r) * NE - r * SW + r * SE
        x13 = (1 - r) * (1 - t) * NBs + (1 - r) * t * NFs + r * (1 - t) * SBs + r * t * SFs
        x23 = -t * WF - (1 - t) * WB + (1 - t) * EB + t * EF
        x123 = [None, None, None]
        for i in range(3):
            x123[i] = - (1 - r) * (1 - t) * NWB[i] \
                      - (1 - r) * t * NWF[i] \
                      + (1 - r) * (1 - t) * NEB[i] \
                      + (1 - r) * t * NEF[i] \
                      - r * (1 - t) * SWB[i] \
                      - r * t * SWF[i] \
                      + r * (1 - t) * SEB[i] \
                      + r * t * SEF[i]
        xs, ys, zs = x1 + x2 + x3 - x12 - x13 - x23 + x123

        x1 = (1 - r) * Nt + r * St
        x2 = (1 - s) * Wt + s * Et
        x3 = -_1_ * B + _1_ * F
        x12 = (1 - r) * (1 - s) * NWt + (1 - r) * s * NEt + r * (1 - s) * SWt + r * s * SEt
        x13 = -(1 - r) * NB + (1 - r) * NF - r * SB + r * SF
        x23 = -(1 - s) * WB + (1 - s) * WF - s * EB + s * EF
        x123 = [None, None, None]
        for i in range(3):
            x123[i] = - (1 - r) * (1 - s) * NWB[i] \
                      + (1 - r) * (1 - s) * NWF[i] \
                      - (1 - r) * s * NEB[i] \
                      + (1 - r) * s * NEF[i] \
                      - r * (1 - s) * SWB[i] \
                      + r * (1 - s) * SWF[i] \
                      - r * s * SEB[i] \
                      + r * s * SEF[i]
        xt, yt, zt = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return ((xr, xs, xt),
                (yr, ys, yt),
                (zr, zs, zt))

    def Jacobian_Xr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        N = self.___NnN___(s, t)[0]
        S = self.___SsS___(s, t)[0]
        Wr = self.___Wr_Wr___(r, t)[0]
        Er = self.___Er_Er___(r, t)[0]
        Br = self.___Br_Br___(r, s)[0]
        Fr = self.___Fr_Fr___(r, s)[0]

        NW = self.___NW___(t)[0]
        NE = self.___NE___(t)[0]
        SW = self.___SW___(t)[0]
        SE = self.___SE___(t)[0]

        NB = self.___NB___(s)[0]
        NF = self.___NF___(s)[0]
        SB = self.___SB___(s)[0]
        SF = self.___SF___(s)[0]
        WBr = self.___WB_r___(r)[0]
        WFr = self.___WF_r___(r)[0]
        EBr = self.___EB_r___(r)[0]
        EFr = self.___EF_r___(r)[0]

        NWB = self.___N_W_B___[0]
        NWF = self.___N_W_F___[0]
        NEB = self.___N_E_B___[0]
        NEF = self.___N_E_F___[0]
        SWB = self.___S_W_B___[0]
        SWF = self.___S_W_F___[0]
        SEB = self.___S_E_B___[0]
        SEF = self.___S_E_F___[0]

        x1 = -_1_ * N + _1_ * S
        x2 = (1 - s) * Wr + s * Er
        x3 = (1 - t) * Br + t * Fr
        x12 = -s * NE - (1 - s) * NW + (1 - s) * SW + s * SE
        x13 = -t * NF - (1 - t) * NB + (1 - t) * SB + t * SF
        x23 = (1 - s) * (1 - t) * WBr + (1 - s) * t * WFr + s * (1 - t) * EBr + s * t * EFr

        x123 = - (1 - s) * (1 - t) * NWB \
               - (1 - s) * t * NWF \
               - s * (1 - t) * NEB \
               - s * t * NEF \
               + (1 - s) * (1 - t) * SWB \
               + (1 - s) * t * SWF \
               + s * (1 - t) * SEB \
               + s * t * SEF
        xr = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return xr

    def Jacobian_Xs(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        Ns = self.___Ns_Ns___(s, t)[0]
        Ss = self.___Ss_Ss___(s, t)[0]

        W = self.___WwW___(r, t)[0]
        E = self.___EeE___(r, t)[0]

        Bs = self.___Bs_Bs___(r, s)[0]
        Fs = self.___Fs_Fs___(r, s)[0]

        NW = self.___NW___(t)[0]
        NE = self.___NE___(t)[0]
        SW = self.___SW___(t)[0]
        SE = self.___SE___(t)[0]

        NBs = self.___NB_s___(s)[0]
        NFs = self.___NF_s___(s)[0]
        SBs = self.___SB_s___(s)[0]
        SFs = self.___SF_s___(s)[0]

        WB = self.___WB___(r)[0]
        WF = self.___WF___(r)[0]
        EB = self.___EB___(r)[0]
        EF = self.___EF___(r)[0]

        NWB = self.___N_W_B___[0]
        NWF = self.___N_W_F___[0]
        NEB = self.___N_E_B___[0]
        NEF = self.___N_E_F___[0]
        SWB = self.___S_W_B___[0]
        SWF = self.___S_W_F___[0]
        SEB = self.___S_E_B___[0]
        SEF = self.___S_E_F___[0]

        x1 = (1 - r) * Ns + r * Ss
        x2 = -_1_ * W + _1_ * E
        x3 = (1 - t) * Bs + t * Fs
        x12 = -(1 - r) * NW + (1 - r) * NE - r * SW + r * SE
        x13 = (1 - r) * (1 - t) * NBs + (1 - r) * t * NFs + r * (1 - t) * SBs + r * t * SFs
        x23 = -t * WF - (1 - t) * WB + (1 - t) * EB + t * EF

        x123 = - (1 - r) * (1 - t) * NWB \
               - (1 - r) * t * NWF \
               + (1 - r) * (1 - t) * NEB \
               + (1 - r) * t * NEF \
               - r * (1 - t) * SWB \
               - r * t * SWF \
               + r * (1 - t) * SEB \
               + r * t * SEF
        xs = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return xs

    def Jacobian_Xt(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))
        Nt = self.___Nt_Nt___(s, t)[0]
        St = self.___St_St___(s, t)[0]

        Wt = self.___Wt_Wt___(r, t)[0]
        Et = self.___Et_Et___(r, t)[0]

        B = self.___BbB___(r, s)[0]
        F = self.___FfF___(r, s)[0]

        NWt = self.___NW_t___(t)[0]
        NEt = self.___NE_t___(t)[0]
        SWt = self.___SW_t___(t)[0]
        SEt = self.___SE_t___(t)[0]

        NB = self.___NB___(s)[0]
        NF = self.___NF___(s)[0]
        SB = self.___SB___(s)[0]
        SF = self.___SF___(s)[0]

        WB = self.___WB___(r)[0]
        WF = self.___WF___(r)[0]
        EB = self.___EB___(r)[0]
        EF = self.___EF___(r)[0]

        NWB = self.___N_W_B___[0]
        NWF = self.___N_W_F___[0]
        NEB = self.___N_E_B___[0]
        NEF = self.___N_E_F___[0]
        SWB = self.___S_W_B___[0]
        SWF = self.___S_W_F___[0]
        SEB = self.___S_E_B___[0]
        SEF = self.___S_E_F___[0]

        x1 = (1 - r) * Nt + r * St
        x2 = (1 - s) * Wt + s * Et
        x3 = -_1_ * B + _1_ * F
        x12 = (1 - r) * (1 - s) * NWt + (1 - r) * s * NEt + r * (1 - s) * SWt + r * s * SEt
        x13 = -(1 - r) * NB + (1 - r) * NF - r * SB + r * SF
        x23 = -(1 - s) * WB + (1 - s) * WF - s * EB + s * EF

        x123 = - (1 - r) * (1 - s) * NWB \
               + (1 - r) * (1 - s) * NWF \
               - (1 - r) * s * NEB \
               + (1 - r) * s * NEF \
               - r * (1 - s) * SWB \
               + r * (1 - s) * SWF \
               - r * s * SEB \
               + r * s * SEF
        xt = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return xt

    def Jacobian_Yr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        N = self.___NnN___(s, t)[1]
        S = self.___SsS___(s, t)[1]

        Wr = self.___Wr_Wr___(r, t)[1]
        Er = self.___Er_Er___(r, t)[1]

        Br = self.___Br_Br___(r, s)[1]
        Fr = self.___Fr_Fr___(r, s)[1]

        NW = self.___NW___(t)[1]
        NE = self.___NE___(t)[1]
        SW = self.___SW___(t)[1]
        SE = self.___SE___(t)[1]

        NB = self.___NB___(s)[1]
        NF = self.___NF___(s)[1]
        SB = self.___SB___(s)[1]
        SF = self.___SF___(s)[1]
        WBr = self.___WB_r___(r)[1]
        WFr = self.___WF_r___(r)[1]
        EBr = self.___EB_r___(r)[1]
        EFr = self.___EF_r___(r)[1]

        NWB = self.___N_W_B___[1]
        NWF = self.___N_W_F___[1]
        NEB = self.___N_E_B___[1]
        NEF = self.___N_E_F___[1]
        SWB = self.___S_W_B___[1]
        SWF = self.___S_W_F___[1]
        SEB = self.___S_E_B___[1]
        SEF = self.___S_E_F___[1]

        x1 = -_1_ * N + _1_ * S
        x2 = (1 - s) * Wr + s * Er
        x3 = (1 - t) * Br + t * Fr
        x12 = -s * NE - (1 - s) * NW + (1 - s) * SW + s * SE
        x13 = -t * NF - (1 - t) * NB + (1 - t) * SB + t * SF
        x23 = (1 - s) * (1 - t) * WBr + (1 - s) * t * WFr + s * (1 - t) * EBr + s * t * EFr

        x123 = - (1 - s) * (1 - t) * NWB \
               - (1 - s) * t * NWF \
               - s * (1 - t) * NEB \
               - s * t * NEF \
               + (1 - s) * (1 - t) * SWB \
               + (1 - s) * t * SWF \
               + s * (1 - t) * SEB \
               + s * t * SEF
        yr = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return yr

    def Jacobian_Ys(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        Ns = self.___Ns_Ns___(s, t)[1]
        Ss = self.___Ss_Ss___(s, t)[1]

        W = self.___WwW___(r, t)[1]
        E = self.___EeE___(r, t)[1]

        Bs = self.___Bs_Bs___(r, s)[1]
        Fs = self.___Fs_Fs___(r, s)[1]

        NW = self.___NW___(t)[1]
        NE = self.___NE___(t)[1]
        SW = self.___SW___(t)[1]
        SE = self.___SE___(t)[1]

        NBs = self.___NB_s___(s)[1]
        NFs = self.___NF_s___(s)[1]
        SBs = self.___SB_s___(s)[1]
        SFs = self.___SF_s___(s)[1]

        WB = self.___WB___(r)[1]
        WF = self.___WF___(r)[1]
        EB = self.___EB___(r)[1]
        EF = self.___EF___(r)[1]

        NWB = self.___N_W_B___[1]
        NWF = self.___N_W_F___[1]
        NEB = self.___N_E_B___[1]
        NEF = self.___N_E_F___[1]
        SWB = self.___S_W_B___[1]
        SWF = self.___S_W_F___[1]
        SEB = self.___S_E_B___[1]
        SEF = self.___S_E_F___[1]

        x1 = (1 - r) * Ns + r * Ss
        x2 = -_1_ * W + _1_ * E
        x3 = (1 - t) * Bs + t * Fs
        x12 = -(1 - r) * NW + (1 - r) * NE - r * SW + r * SE
        x13 = (1 - r) * (1 - t) * NBs + (1 - r) * t * NFs + r * (1 - t) * SBs + r * t * SFs
        x23 = -t * WF - (1 - t) * WB + (1 - t) * EB + t * EF

        x123 = - (1 - r) * (1 - t) * NWB \
               - (1 - r) * t * NWF \
               + (1 - r) * (1 - t) * NEB \
               + (1 - r) * t * NEF \
               - r * (1 - t) * SWB \
               - r * t * SWF \
               + r * (1 - t) * SEB \
               + r * t * SEF
        ys = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return ys

    def Jacobian_Yt(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        Nt = self.___Nt_Nt___(s, t)[1]
        St = self.___St_St___(s, t)[1]

        Wt = self.___Wt_Wt___(r, t)[1]
        Et = self.___Et_Et___(r, t)[1]

        B = self.___BbB___(r, s)[1]
        F = self.___FfF___(r, s)[1]

        NWt = self.___NW_t___(t)[1]
        NEt = self.___NE_t___(t)[1]
        SWt = self.___SW_t___(t)[1]
        SEt = self.___SE_t___(t)[1]

        NB = self.___NB___(s)[1]
        NF = self.___NF___(s)[1]
        SB = self.___SB___(s)[1]
        SF = self.___SF___(s)[1]

        WB = self.___WB___(r)[1]
        WF = self.___WF___(r)[1]
        EB = self.___EB___(r)[1]
        EF = self.___EF___(r)[1]

        NWB = self.___N_W_B___[1]
        NWF = self.___N_W_F___[1]
        NEB = self.___N_E_B___[1]
        NEF = self.___N_E_F___[1]
        SWB = self.___S_W_B___[1]
        SWF = self.___S_W_F___[1]
        SEB = self.___S_E_B___[1]
        SEF = self.___S_E_F___[1]

        x1 = (1 - r) * Nt + r * St
        x2 = (1 - s) * Wt + s * Et
        x3 = -_1_ * B + _1_ * F
        x12 = (1 - r) * (1 - s) * NWt + (1 - r) * s * NEt + r * (1 - s) * SWt + r * s * SEt
        x13 = -(1 - r) * NB + (1 - r) * NF - r * SB + r * SF
        x23 = -(1 - s) * WB + (1 - s) * WF - s * EB + s * EF

        x123 = - (1 - r) * (1 - s) * NWB \
               + (1 - r) * (1 - s) * NWF \
               - (1 - r) * s * NEB \
               + (1 - r) * s * NEF \
               - r * (1 - s) * SWB \
               + r * (1 - s) * SWF \
               - r * s * SEB \
               + r * s * SEF
        yt = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return yt

    def Jacobian_Zr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        N = self.___NnN___(s, t)[2]
        S = self.___SsS___(s, t)[2]

        Wr = self.___Wr_Wr___(r, t)[2]
        Er = self.___Er_Er___(r, t)[2]

        Br = self.___Br_Br___(r, s)[2]
        Fr = self.___Fr_Fr___(r, s)[2]

        NW = self.___NW___(t)[2]
        NE = self.___NE___(t)[2]
        SW = self.___SW___(t)[2]
        SE = self.___SE___(t)[2]

        NB = self.___NB___(s)[2]
        NF = self.___NF___(s)[2]
        SB = self.___SB___(s)[2]
        SF = self.___SF___(s)[2]
        WBr = self.___WB_r___(r)[2]
        WFr = self.___WF_r___(r)[2]
        EBr = self.___EB_r___(r)[2]
        EFr = self.___EF_r___(r)[2]

        NWB = self.___N_W_B___[2]
        NWF = self.___N_W_F___[2]
        NEB = self.___N_E_B___[2]
        NEF = self.___N_E_F___[2]
        SWB = self.___S_W_B___[2]
        SWF = self.___S_W_F___[2]
        SEB = self.___S_E_B___[2]
        SEF = self.___S_E_F___[2]

        x1 = -_1_ * N + _1_ * S
        x2 = (1 - s) * Wr + s * Er
        x3 = (1 - t) * Br + t * Fr
        x12 = -s * NE - (1 - s) * NW + (1 - s) * SW + s * SE
        x13 = -t * NF - (1 - t) * NB + (1 - t) * SB + t * SF
        x23 = (1 - s) * (1 - t) * WBr + (1 - s) * t * WFr + s * (1 - t) * EBr + s * t * EFr

        x123 = - (1 - s) * (1 - t) * NWB \
               - (1 - s) * t * NWF \
               - s * (1 - t) * NEB \
               - s * t * NEF \
               + (1 - s) * (1 - t) * SWB \
               + (1 - s) * t * SWF \
               + s * (1 - t) * SEB \
               + s * t * SEF
        zr = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return zr

    def Jacobian_Zs(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        Ns = self.___Ns_Ns___(s, t)[2]
        Ss = self.___Ss_Ss___(s, t)[2]

        W = self.___WwW___(r, t)[2]
        E = self.___EeE___(r, t)[2]

        Bs = self.___Bs_Bs___(r, s)[2]
        Fs = self.___Fs_Fs___(r, s)[2]

        NW = self.___NW___(t)[2]
        NE = self.___NE___(t)[2]
        SW = self.___SW___(t)[2]
        SE = self.___SE___(t)[2]
        NBs = self.___NB_s___(s)[2]
        NFs = self.___NF_s___(s)[2]
        SBs = self.___SB_s___(s)[2]
        SFs = self.___SF_s___(s)[2]

        WB = self.___WB___(r)[2]
        WF = self.___WF___(r)[2]
        EB = self.___EB___(r)[2]
        EF = self.___EF___(r)[2]

        NWB = self.___N_W_B___[2]
        NWF = self.___N_W_F___[2]
        NEB = self.___N_E_B___[2]
        NEF = self.___N_E_F___[2]
        SWB = self.___S_W_B___[2]
        SWF = self.___S_W_F___[2]
        SEB = self.___S_E_B___[2]
        SEF = self.___S_E_F___[2]

        x1 = (1 - r) * Ns + r * Ss
        x2 = -_1_ * W + _1_ * E
        x3 = (1 - t) * Bs + t * Fs
        x12 = -(1 - r) * NW + (1 - r) * NE - r * SW + r * SE
        x13 = (1 - r) * (1 - t) * NBs + (1 - r) * t * NFs + r * (1 - t) * SBs + r * t * SFs
        x23 = -t * WF - (1 - t) * WB + (1 - t) * EB + t * EF
        x123 = - (1 - r) * (1 - t) * NWB \
               - (1 - r) * t * NWF \
               + (1 - r) * (1 - t) * NEB \
               + (1 - r) * t * NEF \
               - r * (1 - t) * SWB \
               - r * t * SWF \
               + r * (1 - t) * SEB \
               + r * t * SEF
        zs = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return zs

    def Jacobian_Zt(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        _1_ = np.ones(np.shape(r))

        Nt = self.___Nt_Nt___(s, t)[2]
        St = self.___St_St___(s, t)[2]

        Wt = self.___Wt_Wt___(r, t)[2]
        Et = self.___Et_Et___(r, t)[2]

        B = self.___BbB___(r, s)[2]
        F = self.___FfF___(r, s)[2]

        NWt = self.___NW_t___(t)[2]
        NEt = self.___NE_t___(t)[2]
        SWt = self.___SW_t___(t)[2]
        SEt = self.___SE_t___(t)[2]

        NB = self.___NB___(s)[2]
        NF = self.___NF___(s)[2]
        SB = self.___SB___(s)[2]
        SF = self.___SF___(s)[2]

        WB = self.___WB___(r)[2]
        WF = self.___WF___(r)[2]
        EB = self.___EB___(r)[2]
        EF = self.___EF___(r)[2]

        NWB = self.___N_W_B___[2]
        NWF = self.___N_W_F___[2]
        NEB = self.___N_E_B___[2]
        NEF = self.___N_E_F___[2]
        SWB = self.___S_W_B___[2]
        SWF = self.___S_W_F___[2]
        SEB = self.___S_E_B___[2]
        SEF = self.___S_E_F___[2]

        x1 = (1 - r) * Nt + r * St
        x2 = (1 - s) * Wt + s * Et
        x3 = -_1_ * B + _1_ * F
        x12 = (1 - r) * (1 - s) * NWt + (1 - r) * s * NEt + r * (1 - s) * SWt + r * s * SEt
        x13 = -(1 - r) * NB + (1 - r) * NF - r * SB + r * SF
        x23 = -(1 - s) * WB + (1 - s) * WF - s * EB + s * EF

        x123 = - (1 - r) * (1 - s) * NWB \
               + (1 - r) * (1 - s) * NWF \
               - (1 - r) * s * NEB \
               + (1 - r) * s * NEF \
               - r * (1 - s) * SWB \
               + r * (1 - s) * SWF \
               - r * s * SEB \
               + r * s * SEF
        zt = x1 + x2 + x3 - x12 - x13 - x23 + x123

        return zt

    def Jacobian_X_(self, r, s, t):
        return self.Jacobian_Xr(r, s, t), self.Jacobian_Xs(r, s, t), self.Jacobian_Xt(r, s, t)

    def Jacobian_Y_(self, r, s, t):
        return self.Jacobian_Yr(r, s, t), self.Jacobian_Ys(r, s, t), self.Jacobian_Yt(r, s, t)

    def Jacobian_Z_(self, r, s, t):
        return self.Jacobian_Zr(r, s, t), self.Jacobian_Zs(r, s, t), self.Jacobian_Zt(r, s, t)




    def ___NnN___(self, s, t):
        return (self._region_._side_geometries_['N'].X(s, t),
                self._region_._side_geometries_['N'].Y(s, t),
                self._region_._side_geometries_['N'].Z(s, t))

    def ___Ns_Ns___(self, s, t):
        return (self._region_._side_geometries_['N'].Xp(s, t),
                self._region_._side_geometries_['N'].Yp(s, t),
                self._region_._side_geometries_['N'].Zp(s, t))

    def ___Nt_Nt___(self, s, t):
        return (self._region_._side_geometries_['N'].Xq(s, t),
                self._region_._side_geometries_['N'].Yq(s, t),
                self._region_._side_geometries_['N'].Zq(s, t))

    def ___SsS___(self, s, t):
        return (self._region_._side_geometries_['S'].X(s, t),
                self._region_._side_geometries_['S'].Y(s, t),
                self._region_._side_geometries_['S'].Z(s, t))

    def ___Ss_Ss___(self, s, t):
        return (self._region_._side_geometries_['S'].Xp(s, t),
                self._region_._side_geometries_['S'].Yp(s, t),
                self._region_._side_geometries_['S'].Zp(s, t))

    def ___St_St___(self, s, t):
        return (self._region_._side_geometries_['S'].Xq(s, t),
                self._region_._side_geometries_['S'].Yq(s, t),
                self._region_._side_geometries_['S'].Zq(s, t))

    def ___WwW___(self, r, t):
        return (self._region_._side_geometries_['W'].X(r, t),
                self._region_._side_geometries_['W'].Y(r, t),
                self._region_._side_geometries_['W'].Z(r, t))

    def ___Wr_Wr___(self, r, t):
        return (self._region_._side_geometries_['W'].Xp(r, t),
                self._region_._side_geometries_['W'].Yp(r, t),
                self._region_._side_geometries_['W'].Zp(r, t))

    def ___Wt_Wt___(self, r, t):
        return (self._region_._side_geometries_['W'].Xq(r, t),
                self._region_._side_geometries_['W'].Yq(r, t),
                self._region_._side_geometries_['W'].Zq(r, t))

    def ___EeE___(self, r, t):
        return (self._region_._side_geometries_['E'].X(r, t),
                self._region_._side_geometries_['E'].Y(r, t),
                self._region_._side_geometries_['E'].Z(r, t))

    def ___Er_Er___(self, r, t):
        return (self._region_._side_geometries_['E'].Xp(r, t),
                self._region_._side_geometries_['E'].Yp(r, t),
                self._region_._side_geometries_['E'].Zp(r, t))

    def ___Et_Et___(self, r, t):
        return (self._region_._side_geometries_['E'].Xq(r, t),
                self._region_._side_geometries_['E'].Yq(r, t),
                self._region_._side_geometries_['E'].Zq(r, t))

    def ___BbB___(self, r, s):
        return (self._region_._side_geometries_['B'].X(r, s),
                self._region_._side_geometries_['B'].Y(r, s),
                self._region_._side_geometries_['B'].Z(r, s))

    def ___Br_Br___(self, r, s):
        return (self._region_._side_geometries_['B'].Xp(r, s),
                self._region_._side_geometries_['B'].Yp(r, s),
                self._region_._side_geometries_['B'].Zp(r, s))

    def ___Bs_Bs___(self, r, s):
        return (self._region_._side_geometries_['B'].Xq(r, s),
                self._region_._side_geometries_['B'].Yq(r, s),
                self._region_._side_geometries_['B'].Zq(r, s))

    def ___FfF___(self, r, s):
        return (self._region_._side_geometries_['F'].X(r, s),
                self._region_._side_geometries_['F'].Y(r, s),
                self._region_._side_geometries_['F'].Z(r, s))

    def ___Fr_Fr___(self, r, s):
        return (self._region_._side_geometries_['F'].Xp(r, s),
                self._region_._side_geometries_['F'].Yp(r, s),
                self._region_._side_geometries_['F'].Zp(r, s))

    def ___Fs_Fs___(self, r, s):
        return (self._region_._side_geometries_['F'].Xq(r, s),
                self._region_._side_geometries_['F'].Yq(r, s),
                self._region_._side_geometries_['F'].Zq(r, s))

    def ___NW___(self, t):
        return (self._region_._side_geometries_['N'].X(0, t),
                self._region_._side_geometries_['N'].Y(0, t),
                self._region_._side_geometries_['N'].Z(0, t))

    def ___NW_t___(self, t):
        return (self._region_._side_geometries_['N'].Xq(0, t),
                self._region_._side_geometries_['N'].Yq(0, t),
                self._region_._side_geometries_['N'].Zq(0, t))

    def ___NE___(self, t):
        return (self._region_._side_geometries_['N'].X(1, t),
                self._region_._side_geometries_['N'].Y(1, t),
                self._region_._side_geometries_['N'].Z(1, t))

    def ___NE_t___(self, t):
        return (self._region_._side_geometries_['N'].Xq(1, t),
                self._region_._side_geometries_['N'].Yq(1, t),
                self._region_._side_geometries_['N'].Zq(1, t))

    def ___NB___(self, s):
        return (self._region_._side_geometries_['N'].X(s, 0),
                self._region_._side_geometries_['N'].Y(s, 0),
                self._region_._side_geometries_['N'].Z(s, 0))

    def ___NB_s___(self, s):
        return (self._region_._side_geometries_['N'].Xp(s, 0),
                self._region_._side_geometries_['N'].Yp(s, 0),
                self._region_._side_geometries_['N'].Zp(s, 0))

    def ___NF___(self, s):
        return (self._region_._side_geometries_['N'].X(s, 1),
                self._region_._side_geometries_['N'].Y(s, 1),
                self._region_._side_geometries_['N'].Z(s, 1))

    def ___NF_s___(self, s):
        return (self._region_._side_geometries_['N'].Xp(s, 1),
                self._region_._side_geometries_['N'].Yp(s, 1),
                self._region_._side_geometries_['N'].Zp(s, 1))

    def ___SW___(self, t):
        return (self._region_._side_geometries_['S'].X(0, t),
                self._region_._side_geometries_['S'].Y(0, t),
                self._region_._side_geometries_['S'].Z(0, t))

    def ___SW_t___(self, t):
        return (self._region_._side_geometries_['S'].Xq(0, t),
                self._region_._side_geometries_['S'].Yq(0, t),
                self._region_._side_geometries_['S'].Zq(0, t))

    def ___SE___(self, t):
        return (self._region_._side_geometries_['S'].X(1, t),
                self._region_._side_geometries_['S'].Y(1, t),
                self._region_._side_geometries_['S'].Z(1, t))

    def ___SE_t___(self, t):
        return (self._region_._side_geometries_['S'].Xq(1, t),
                self._region_._side_geometries_['S'].Yq(1, t),
                self._region_._side_geometries_['S'].Zq(1, t))

    def ___SB___(self, s):
        return (self._region_._side_geometries_['S'].X(s, 0),
                self._region_._side_geometries_['S'].Y(s, 0),
                self._region_._side_geometries_['S'].Z(s, 0))

    def ___SB_s___(self, s):
        return (self._region_._side_geometries_['S'].Xp(s, 0),
                self._region_._side_geometries_['S'].Yp(s, 0),
                self._region_._side_geometries_['S'].Zp(s, 0))

    def ___SF___(self, s):
        return (self._region_._side_geometries_['S'].X(s, 1),
                self._region_._side_geometries_['S'].Y(s, 1),
                self._region_._side_geometries_['S'].Z(s, 1))

    def ___SF_s___(self, s):
        return (self._region_._side_geometries_['S'].Xp(s, 1),
                self._region_._side_geometries_['S'].Yp(s, 1),
                self._region_._side_geometries_['S'].Zp(s, 1))

    def ___WN___(self, t):
        return (self._region_._side_geometries_['W'].X(0, t),
                self._region_._side_geometries_['W'].Y(0, t),
                self._region_._side_geometries_['W'].Z(0, t))

    def ___WS___(self, t):
        return (self._region_._side_geometries_['W'].X(1, t),
                self._region_._side_geometries_['W'].Y(1, t),
                self._region_._side_geometries_['W'].Z(1, t))

    def ___WB___(self, r):
        return (self._region_._side_geometries_['W'].X(r, 0),
                self._region_._side_geometries_['W'].Y(r, 0),
                self._region_._side_geometries_['W'].Z(r, 0))

    def ___WB_r___(self, r):
        return (self._region_._side_geometries_['W'].Xp(r, 0),
                self._region_._side_geometries_['W'].Yp(r, 0),
                self._region_._side_geometries_['W'].Zp(r, 0))

    def ___WF___(self, r):
        return (self._region_._side_geometries_['W'].X(r, 1),
                self._region_._side_geometries_['W'].Y(r, 1),
                self._region_._side_geometries_['W'].Z(r, 1))

    def ___WF_r___(self, r):
        return (self._region_._side_geometries_['W'].Xp(r, 1),
                self._region_._side_geometries_['W'].Yp(r, 1),
                self._region_._side_geometries_['W'].Zp(r, 1))

    def ___EN___(self, t):
        return (self._region_._side_geometries_['E'].X(0, t),
                self._region_._side_geometries_['E'].Y(0, t),
                self._region_._side_geometries_['E'].Z(0, t))

    def ___ES___(self, t):
        return (self._region_._side_geometries_['E'].X(1, t),
                self._region_._side_geometries_['E'].Y(1, t),
                self._region_._side_geometries_['E'].Z(1, t))

    def ___EB___(self, r):
        return (self._region_._side_geometries_['E'].X(r, 0),
                self._region_._side_geometries_['E'].Y(r, 0),
                self._region_._side_geometries_['E'].Z(r, 0))

    def ___EB_r___(self, r):
        return (self._region_._side_geometries_['E'].Xp(r, 0),
                self._region_._side_geometries_['E'].Yp(r, 0),
                self._region_._side_geometries_['E'].Zp(r, 0))

    def ___EF___(self, r):
        return (self._region_._side_geometries_['E'].X(r, 1),
                self._region_._side_geometries_['E'].Y(r, 1),
                self._region_._side_geometries_['E'].Z(r, 1))

    def ___EF_r___(self, r):
        return (self._region_._side_geometries_['E'].Xp(r, 1),
                self._region_._side_geometries_['E'].Yp(r, 1),
                self._region_._side_geometries_['E'].Zp(r, 1))

    def ___BN___(self, s):
        return (self._region_._side_geometries_['B'].X(0, s),
                self._region_._side_geometries_['B'].Y(0, s),
                self._region_._side_geometries_['B'].Z(0, s))

    def ___BS___(self, s):
        return (self._region_._side_geometries_['B'].X(1, s),
                self._region_._side_geometries_['B'].Y(1, s),
                self._region_._side_geometries_['B'].Z(1, s))

    def ___BW___(self, r):
        return (self._region_._side_geometries_['B'].X(r, 0),
                self._region_._side_geometries_['B'].Y(r, 0),
                self._region_._side_geometries_['B'].Z(r, 0))

    def ___BE___(self, r):
        return (self._region_._side_geometries_['B'].X(r, 1),
                self._region_._side_geometries_['B'].Y(r, 1),
                self._region_._side_geometries_['B'].Z(r, 1))

    def ___FN___(self, s):
        return (self._region_._side_geometries_['F'].X(0, s),
                self._region_._side_geometries_['F'].Y(0, s),
                self._region_._side_geometries_['F'].Z(0, s))

    def ___FS___(self, s):
        return (self._region_._side_geometries_['F'].X(1, s),
                self._region_._side_geometries_['F'].Y(1, s),
                self._region_._side_geometries_['F'].Z(1, s))

    def ___FW___(self, r):
        return (self._region_._side_geometries_['F'].X(r, 0),
                self._region_._side_geometries_['F'].Y(r, 0),
                self._region_._side_geometries_['F'].Z(r, 0))

    def ___FE___(self, r):
        return (self._region_._side_geometries_['F'].X(r, 1),
                self._region_._side_geometries_['F'].Y(r, 1),
                self._region_._side_geometries_['F'].Z(r, 1))

    @property
    def ___N_W_B___(self):
        return (self._region_._side_geometries_['N'].X(0, 0),
                self._region_._side_geometries_['N'].Y(0, 0),
                self._region_._side_geometries_['N'].Z(0, 0))

    @property
    def ___N_W_F___(self):
        return (self._region_._side_geometries_['N'].X(0, 1),
                self._region_._side_geometries_['N'].Y(0, 1),
                self._region_._side_geometries_['N'].Z(0, 1))

    @property
    def ___N_E_B___(self):
        return (self._region_._side_geometries_['N'].X(1, 0),
                self._region_._side_geometries_['N'].Y(1, 0),
                self._region_._side_geometries_['N'].Z(1, 0))

    @property
    def ___N_E_F___(self):
        return (self._region_._side_geometries_['N'].X(1, 1),
                self._region_._side_geometries_['N'].Y(1, 1),
                self._region_._side_geometries_['N'].Z(1, 1))

    @property
    def ___S_W_B___(self):
        return (self._region_._side_geometries_['S'].X(0, 0),
                self._region_._side_geometries_['S'].Y(0, 0),
                self._region_._side_geometries_['S'].Z(0, 0))

    @property
    def ___S_W_F___(self):
        return (self._region_._side_geometries_['S'].X(0, 1),
                self._region_._side_geometries_['S'].Y(0, 1),
                self._region_._side_geometries_['S'].Z(0, 1))

    @property
    def ___S_E_B___(self):
        return (self._region_._side_geometries_['S'].X(1, 0),
                self._region_._side_geometries_['S'].Y(1, 0),
                self._region_._side_geometries_['S'].Z(1, 0))

    @property
    def ___S_E_F___(self):
        return (self._region_._side_geometries_['S'].X(1, 1),
                self._region_._side_geometries_['S'].Y(1, 1),
                self._region_._side_geometries_['S'].Z(1, 1))