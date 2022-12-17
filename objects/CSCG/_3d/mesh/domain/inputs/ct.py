# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/30/2022 2:35 PM
"""
from objects.CSCG._3d.mesh.domain.inputs.base import _3dDomainInputBase
from components.decorators.classproperty.main import classproperty

import random
from root.config.main import RANK, MASTER_RANK, COMM


class CurvilinearTestMesh(_3dDomainInputBase):
    """A curvilinear mesh for test purpose."""

    def __init__(self, c=0.1):
        """"""
        self._c_ = c
        super().__init__(domain_name='CurvilinearTest')

        x0, x1 = [0, 1]
        y0, y1 = [0, 1]
        z0, z1 = [0, 1]

        assert x1 > x0
        assert y1 > y0
        assert z1 > z0
        self.region_corner_coordinates = {'R:R': ((x0, y0, z0), (x1, y0, z0), (x0, y1, z0), (x1, y1, z0),
                                                  (x0, y0, z1), (x1, y0, z1), (x0, y1, z1), (x1, y1, z1))}
        self.region_side_types = {'R:R-S': ('free',)}
        self.boundary_region_sides = {'South': "R:R-S", 'North': "R:R-N",
                                      'West': "R:R-W", 'East': "R:R-E",
                                      'Back': "R:R-B", 'Front': "R:R-F"}
        self.region_interpolators = 'ct'
        self.region_type_wr2_metric = {'R:R': 'chaotic'}
        self.internal_parameters = ['c', ]  # has to be defined after the super().__init__

        self._freeze_self_()

    @property
    def c(self):
        return self._c_

    @classproperty
    def statistic(cls):
        return {
            'periodic': False,
            'region num': 1,
            'mesh boundary num': 6,  # the amount of mesh boundaries (instead of domain boundaries)
        }

    @classproperty
    def random_parameters(cls):
        if RANK == MASTER_RANK:
            rp = {
                'c': random.randint(0, 3) * random.random() / 10,
            }
        else:
            rp = None

        return COMM.bcast(rp, root=MASTER_RANK)
