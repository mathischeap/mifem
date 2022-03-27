
from objects.CSCG._2d.mesh.domain.inputs.base import DomainInputBase
import numpy as np
from screws.decorators.classproperty.main import classproperty

class CircleHolePlate2(DomainInputBase):
    """ """
    def __init__(self, hx=2, hy=None, r=0.5):
        """
        The DomainInput describes such a domain:

              ________hy_________
             |                   |
             |        ___        |hx/2
             |       /   \       |
           hx|      |  .r |      |-----> y
             |       \___/       |
             |                   |
             |___________________|
                 hy/2  |
                       |
                       |
                       v x

        The domain is divided into following regions:
                       U
              ___________________
             |       |R_U|       |
             | R_UL  |___|  R_UR |
             |-------/   \-------|
           L | R_L  |  .  | R_R  | R
             |-------\___/-------|
             | R_DL  |R_D|  R_DR |
             |_______|___|_______|
                       D

        The center of the circle hole is at (0, 0).

        Parameters
        ----------

        """
        # ____ parse inputs ____________________________________________________________
        if hy is None: hy = hx
        # ____ checks __________________________________________________________________
        assert hx > 0 and hy > 0, " <HolePlate> : hx={}, hy={} illegal.".format(hx, hy)
        assert r < np.min((hx, hy)), " <HolePlate> : r={} too large.".format(r)
        # _____________ standard inputs ________________________________________________
        super().__init__(domain_name='CircleHolePlate2')
        sr = np.sqrt(2) * r / 2
        self.region_corner_coordinates = {
            'R:R_UL': ((-hx / 2, -hy / 2), (-sr, -hy / 2), (-hx / 2, -sr), (-sr, -sr)),
            'R:R_L': ((-sr, -hy / 2), (sr, -hy / 2), (-sr, -sr), (sr, -sr)),
            'R:R_DL': ((sr, -hy / 2), (hx / 2, -hy / 2), (sr, -sr), (hx / 2, -sr)),
            'R:R_U': ((-hx / 2, -sr), (-sr, -sr), (-hx / 2, sr), (-sr, sr)),
            'R:R_D': ((sr, -sr), (hx / 2, -sr), (sr, sr), (hx / 2, sr)),
            'R:R_UR': ((-hx / 2, sr), (-sr, sr), (-hx / 2, hy / 2), (-sr, hy / 2)),
            'R:R_R': ((-sr, sr), (sr, sr), (-sr, hy / 2), (sr, hy / 2)),
            'R:R_DR': ((sr, sr), (hx / 2, sr), (sr, hy / 2), (hx / 2, hy / 2))}
        self.region_edge_types = {'R:R_L-R': ('aacw', (0, 0)),
                                  'R:R_U-D': ('acw', (0, 0)),
                                  'R:R_D-U': ('aacw', (0, 0)),
                                  'R:R_R-L': ('acw', (0, 0)), }
        self.boundary_region_edges = {'Upper': ("R:R_UL-U", 'R:R_U-U', "R:R_UR-U"),
                                      'Down': ("R:R_DL-D", 'R:R_D-D', "R:R_DR-D"),
                                      'Left': ("R:R_UL-L", 'R:R_L-L', "R:R_DL-L"),
                                      'Right': ("R:R_UR-R", 'R:R_R-R', "R:R_DR-R"),
                                      'Internal': ("R:R_L-R", "R:R_U-D", "R:R_D-U", 'R:R_R-L')}
        self.region_interpolators = 'transfinite'

        self.region_type_wr2_metric = 'transfinite'
        self.internal_parameters = list()


    @classproperty
    def statistic(cls):
        return {'periodic': False,
                'region num': 8,
                'mesh boundary num': 5, # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        return {}