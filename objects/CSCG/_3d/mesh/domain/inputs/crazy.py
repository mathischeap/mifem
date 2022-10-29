




from objects.CSCG._3d.mesh.domain.inputs.base import _3dDomainInputBase
import numpy as np
from screws.decorators.classproperty.main import classproperty



import random
from root.config.main import RANK, MASTER_RANK, COMM

class Crazy(_3dDomainInputBase):
    """A "crazy" 3d rectangular domain's input class whose inside is distorted with the "crazy" mapping."""

    def __init__(self, c=0, bounds=((0, 1), (0, 1), (0, 1)), domain_name="Crazy"):
        assert np.shape(bounds)[0] == 3, " <CrazyDomain> : bounds dimension is wrong."
        for i in range(3):
            assert np.shape(bounds[i]) == (2,) and bounds[i][1] > bounds[i][0], \
                " <CrazyDomain> : bounds[{}]=={} is wrong.".format(i, bounds[i])
        self._bounds_ = bounds
        self._c_ = c
        super().__init__(domain_name=domain_name)

        x0, x1 = bounds[0]
        y0, y1 = bounds[1]
        z0, z1 = bounds[2]
        assert x1 > x0
        assert y1 > y0
        assert z1 > z0
        self.region_corner_coordinates = {'R:R': ((x0, y0, z0), (x1, y0, z0), (x0, y1, z0), (x1, y1, z0),
                                                  (x0, y0, z1), (x1, y0, z1), (x0, y1, z1), (x1, y1, z1))}
        self.region_side_types = {'R:R-S': ('free',)}
        self.boundary_region_sides = {'South': "R:R-S", 'North': "R:R-N",
                                      'West': "R:R-W", 'East': "R:R-E",
                                      'Back': "R:R-B", 'Front': "R:R-F"}
        self.region_interpolators = 'crazy'
        self.region_type_wr2_metric = {'R:R': 'crazy'}
        self.internal_parameters = ['c', ]  # has to be defined after the super().__init__

    @property
    def bounds(self):
        return self._bounds_

    @property
    def c(self):
        return self._c_

    @classproperty
    def statistic(cls):
        return {'periodic': False,
                'region num': 1,
                'mesh boundary num': 6, # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        if RANK == MASTER_RANK:
            rp = {'c': random.randint(0,3) * random.random() / 10,
                  'bounds': [(-random.random(), random.random()+0.5),
                             (-random.random(), random.random()+0.5),
                             (-random.random(), random.random()+0.5)]
                   }
        else:
            rp = None

        return COMM.bcast(rp, root=MASTER_RANK)