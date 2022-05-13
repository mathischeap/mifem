# -*- coding: utf-8 -*-

from objects.CSCG._2d.mesh.domain.inputs.base import DomainInputBase
import numpy as np
from screws.decorators.classproperty.main import classproperty

class CylinderInChannel(DomainInputBase):
    """ Just like the class name say, this is a cylinder in channel domain."""
    def __init__(self, h=1.5, r=0.5, li=0.75, lo=2.25):
        """
        The DomainInput describes such a domain:
                            y
                            ^
                            |
         ____li_____._______|________.___________lo___________
        |           .       |        .                        |
        |           .     __|__      .                        |
        |           .    /  |  \     .                        |
        |h          .h  |   .-r-|----.------------------------|-----> x
        |           .    \_____/     .                        |
        |           .                .                        |
        |____li_____._______h________.___________lo___________|

        The domain is divided into following regions:
         ___________._______________.________________________
        |           .\      Ru     /.                        |
        |           . \   _____   / .                        |
        |           .  \ /     \ /  .                        |
        |     Ri    .Rl |   .-r-| Rr.           Ro           |
        |           .  / \_____/ \  .                        |
        |           . /     Rd    \ .                        |
        |___________./_____________\.________________________|


        Parameters
        ----------
        h : float
            The height. It must > 0.
        r : float
            The radius of the cylinder. It must > 0.
        li :
            The Inlet Length. It must > 0.
        lo :
            The Outlet Length. It must > 0.

        """
        # ____ do some checks first ____________________________________________________
        assert isinstance(h, (int, float)) and h > 0, " <CylinderInChannel> : need h > 0."
        assert isinstance(r, (int, float)) and r > 0, " <CylinderInChannel> : need r > 0."
        assert isinstance(li, (int, float)) and r > 0, " <CylinderInChannel> : need li > 0."
        assert isinstance(lo, (int, float)) and r > 0, " <CylinderInChannel> : need lo > 0."
        assert r < h / 2, " <CylinderInChannel> : need r < h/2."
        # ____ personal properties _____________________________________________________
        self._h_ = h
        self._r_ = r
        self._li_ = li
        self._lo_ = lo
        # ______ initialiing parent ____________________________________________________
        hsqr_r = 0.5 * np.sqrt(r)
        super().__init__(domain_name='CylinderInChannel')
        self.region_corner_coordinates = {
            'R:Ru': ((-hsqr_r, hsqr_r), (hsqr_r, hsqr_r), (-h / 2, h / 2), (h / 2, h / 2)),
            'R:Rr': ((hsqr_r, hsqr_r), (hsqr_r, -hsqr_r), (h / 2, h / 2), (h / 2, -h / 2)),
            'R:Rd': ((hsqr_r, -hsqr_r), (-hsqr_r, -hsqr_r), (h / 2, -h / 2), (-h / 2, -h / 2)),
            'R:Rl': ((-hsqr_r, -hsqr_r), (-hsqr_r, hsqr_r), (-h / 2, -h / 2), (-h / 2, h / 2)),
            'R:Ri': ((-h / 2, -h / 2), (-h / 2, h / 2), (-h / 2 - li, -h / 2), (-h / 2 - li, h / 2)),
            'R:Ro': ((h / 2, h / 2), (h / 2, -h / 2), (h / 2 + lo, h / 2), (h / 2 + lo, -h / 2))}
        self.region_edge_types = {'R:Ru-L': ('acw', (0, 0)),
                                  'R:Rl-L': ('acw', (0, 0)),
                                  'R:Rd-L': ('acw', (0, 0)),
                                  'R:Rr-L': ('acw', (0, 0)), }
        self.boundary_region_edges = {'Upper': ('R:Ri-D', 'R:Ru-R', 'R:Ro-U'),
                                      'Down': ('R:Ri-U', 'R:Rd-R', 'R:Ro-D'),
                                      'Left': 'R:Ri-R',
                                      'Right': 'R:Ro-R',
                                      'Internal': ('R:Ru-L', 'R:Rl-L', 'R:Rd-L', 'R:Rr-L')}
        self.region_interpolators = 'transfinite'

        self.region_type_wr2_metric = 'transfinite'
        # self.region_sequence = ('R:R',)
        # self._internal_parameters_ = None

        # ------------------------------------------------------------------------------

    @property
    def h(self):
        return self._h_

    @property
    def r(self):
        return self._r_

    @property
    def li(self):
        return self._li_

    @property
    def lo(self):
        return self._lo_




    @classproperty
    def statistic(cls):
        return {'periodic': False,
                'region num': 6,
                'mesh boundary num': 5, # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        return {}