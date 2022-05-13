# -*- coding: utf-8 -*-
from objects.CSCG._2d.mesh.domain.inputs.base import DomainInputBase
import numpy as np
from screws.decorators.classproperty.main import classproperty



class CircleHolePlate1(DomainInputBase):
    """ """
    def __init__(self, hx=1.5, hy=None, r=0.5):
        """
        The DomainInput describes such a domain:
                y
                ^
           hx/2 |  hx/2
         _______|________
        |       |        |
        |     __|__      |hy/2
        |    /  |  \     |
        |hy |   .-r-|----|---> x
        |    \_____/     |
        |                |hy/2
        |_______hx_______|

        The domain is divided into following regions:
         _______________
        |\      Ru     /|
        | \   _____   / |
        |  \ /     \ /  |
        |Rl |   .-r-| Rr|
        |  / \_____/ \  |
        | /     Rd    \ |
        |/_____________\|


        Parameters
        ----------

        """
        # ____ parse inputs ____________________________________________________________
        if hy is None: hy = hx
        # ____ checks __________________________________________________________________
        assert hx > 0 and hy > 0, " <HolePlate> : hx={}, hy={} illegal.".format(hx, hy)
        assert r < np.min((hx, hy)), " <HolePlate> : r={} too large.".format(r)
        # ______ initialiing parent ____________________________________________________
        hsqr_r = 0.5 * np.sqrt(r)
        super().__init__(domain_name='CircleHolePlate1')
        self.region_corner_coordinates = {
            'R:Ru': ((-hsqr_r, hsqr_r), (hsqr_r, hsqr_r), (-hx / 2, hy / 2), (hx / 2, hy / 2)),
            'R:Rr': ((hsqr_r, hsqr_r), (hsqr_r, -hsqr_r), (hx / 2, hy / 2), (hx / 2, -hy / 2)),
            'R:Rd': ((hsqr_r, -hsqr_r), (-hsqr_r, -hsqr_r), (hx / 2, -hy / 2), (-hx / 2, -hy / 2)),
            'R:Rl': ((-hsqr_r, -hsqr_r), (-hsqr_r, hsqr_r), (-hx / 2, -hy / 2), (-hx / 2, hy / 2))}
        self.region_edge_types = {'R:Ru-L': ('acw', (0, 0)),
                                  'R:Rl-L': ('acw', (0, 0)),
                                  'R:Rd-L': ('acw', (0, 0)),
                                  'R:Rr-L': ('acw', (0, 0)), }
        self.boundary_region_edges = {'Upper': 'R:Ru-R',
                                      'Down': 'R:Rd-R',
                                      'Left': 'R:Rl-R',
                                      'Right': 'R:Rr-R',
                                      'Internal': ('R:Ru-L', 'R:Rl-L', 'R:Rd-L', 'R:Rr-L')}
        self.region_interpolators = 'transfinite'

        self.region_type_wr2_metric = 'transfinite'
        self.internal_parameters = list()
        # ------------------------------------------------------------------------------


    @classproperty
    def statistic(cls):
        return {'periodic': False,
                'region num': 4,
                'mesh boundary num': 5, # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        return {}