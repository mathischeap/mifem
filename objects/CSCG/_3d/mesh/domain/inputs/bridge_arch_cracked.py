



from objects.CSCG._3d.mesh.domain.inputs.base import _3dDomainInputBase
import numpy as np
from screws.decorators.classproperty.main import classproperty



class BridgeArchCracked(_3dDomainInputBase):
    def __init__(self, domain_name="BridgeArch",
        l=4, w=1, h=1.5, r=5 / 2, d=0.25,
        region_interpolators='bridge_arch_cracked'):
        """
        Parameters
        ----------
        r : float
            Arch radius.
        d : float
            Uncracked depth.

        """
        self._l_, self._w_, self._h_, self._r_ = l, w, h, r
        assert r ** 2 - l ** 2 / 4 > 0, " <BridgeArch> : parameters do not fit an arch."
        alpha = np.sqrt(r ** 2 - l ** 2 / 4)
        beta = h + alpha - r
        assert beta > 0, " <BridgeArch> : bridge has zero thickness."
        self._center_ = (h + alpha, l / 2)
        self._alpha_ = alpha
        self._beta_ = beta
        self._d_ = d
        assert d < beta, " <BridgeArchCracked> : Bridge already breaks into two."
        A = (0, 0, 0)
        B = (h, 0, 0)
        C = (0, l / 2, 0)
        Dl = (beta, l / 2, 0)
        Dr = (beta, l / 2, 0)
        E = (0, 0, w)
        F = (h, 0, w)
        G = (0, l / 2, w)
        Hl = (beta, l / 2, w)
        Hr = (beta, l / 2, w)
        I = (0, l, 0)
        J = (h, l, 0)
        K = (0, l, w)
        L = (h, l, w)
        rd = beta - d
        M = (rd * beta / h, 0, 0)
        N = (rd * beta / h, 0, w)
        O = (rd, l / 2, 0)
        P = (rd, l / 2, w)
        Q = (rd * beta / h, l, 0)
        R = (rd * beta / h, l, w)
        super().__init__(domain_name=domain_name)
        self.region_corner_coordinates = {'R:R_left_up': (A, M, C, O, E, N, G, P),
                                          'R:R_left_down': (M, B, O, Dl, N, F, P, Hl),
                                          'R:R_right_up': (C, O, I, Q, G, P, K, R),
                                          'R:R_right_down': (O, Dr, Q, J, P, Hr, R, L)}
        self.region_side_types = {} # all plane
        self.boundary_region_sides = {
            'Left_Floor': ("R:R_left_up-N",),
            'Right_Floor': ('R:R_right_up-N',),
            'Back': ("R:R_left_up-B", "R:R_left_down-B", 'R:R_right_up-B', 'R:R_right_down-B',),
            'Front': ("R:R_left_up-F", "R:R_left_down-F", 'R:R_right_up-F', 'R:R_right_down-F',),
            'Bottom': ("R:R_left_down-S", 'R:R_right_down-S',),
            'Left_Wall': ('R:R_left_up-W', 'R:R_left_down-W',),
            'Right_Wall': ('R:R_right_up-E', 'R:R_right_down-E',),
            'Crack': ('R:R_left_down-E', 'R:R_right_down-W')
        }
        self.region_interpolators = region_interpolators

    @property
    def l(self):
        return self._l_

    @property
    def w(self):
        return self._w_

    @property
    def h(self):
        return self._h_

    @property
    def r(self):
        """Arch radius."""
        return self._r_

    @property
    def d(self):
        """Crack depth."""
        return self._d_

    @property
    def center(self):
        return self._center_

    @property
    def beta(self):
        return self._beta_


    @classproperty
    def statistic(cls):
        return {'periodic': False,
                'region num': 4,
                'mesh boundary num': 8, # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        return dict()