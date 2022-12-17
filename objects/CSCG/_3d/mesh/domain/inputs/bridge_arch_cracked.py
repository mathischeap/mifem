# -*- coding: utf-8 -*-
from objects.CSCG._3d.mesh.domain.inputs.base import _3dDomainInputBase
import numpy as np
from components.decorators.classproperty.main import classproperty



class BridgeArchCracked(_3dDomainInputBase):
    def __init__(
            self, domain_name="BridgeArch",
            length=4, w=1, h=1.5, r=5 / 2, d=0.25,
            region_interpolators='bridge_arch_cracked'
    ):
        """
        Parameters
        ----------
        r : float
            Arch radius.
        d : float
            Uncracked depth.

        """
        self._length_, self._w_, self._h_, self._r_ = length, w, h, r
        assert r ** 2 - length ** 2 / 4 > 0, " <BridgeArch> : parameters do not fit an arch."
        alpha = np.sqrt(r ** 2 - length ** 2 / 4)
        beta = h + alpha - r
        assert beta > 0, " <BridgeArch> : bridge has zero thickness."
        self._center_ = (h + alpha, length / 2)
        self._alpha_ = alpha
        self._beta_ = beta
        self._d_ = d
        assert d < beta, " <BridgeArchCracked> : Bridge already breaks into two."
        A = (0, 0, 0)
        B = (h, 0, 0)
        C = (0, length / 2, 0)
        Dl = (beta, length / 2, 0)
        Dr = (beta, length / 2, 0)
        E = (0, 0, w)
        F = (h, 0, w)
        G = (0, length / 2, w)
        Hl = (beta, length / 2, w)
        Hr = (beta, length / 2, w)
        _I = (0, length, 0)
        J = (h, length, 0)
        K = (0, length, w)
        L = (h, length, w)
        rd = beta - d
        M = (rd * beta / h, 0, 0)
        N = (rd * beta / h, 0, w)
        _O = (rd, length / 2, 0)
        P = (rd, length / 2, w)
        Q = (rd * beta / h, length, 0)
        R = (rd * beta / h, length, w)
        super().__init__(domain_name=domain_name)
        self.region_corner_coordinates = {'R:R_left_up': (A, M, C, _O, E, N, G, P),
                                          'R:R_left_down': (M, B, _O, Dl, N, F, P, Hl),
                                          'R:R_right_up': (C, _O, _I, Q, G, P, K, R),
                                          'R:R_right_down': (_O, Dr, Q, J, P, Hr, R, L)}
        self.region_side_types = dict()  # all plane
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
    def L(self):
        return self._length_

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
                'mesh boundary num': 8,  # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        return dict()
