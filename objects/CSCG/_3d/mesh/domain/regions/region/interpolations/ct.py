# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/30/2022 2:56 PM
"""
import numpy as np
from objects.CSCG._3d.mesh.domain.regions.region.interpolations.base import InterpolationBase


class CurvilinearTest(InterpolationBase):
    """"""

    def __init__(self, region):
        """"""
        super(CurvilinearTest, self).__init__(region)
        self._c_ = region._domain_input_.c
        self._freeze_self_()

    def ___PrCT3CSCG_curvature____(self, x, y):
        """"""
        return np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y) * self._c_

    def mapping_X(self, r, s, t):
        return r + self.___PrCT3CSCG_curvature____(s, t)

    def mapping_Y(self, r, s, t):
        return s + self.___PrCT3CSCG_curvature____(r, t)

    def mapping_Z(self, r, s, t):
        return t - self.___PrCT3CSCG_curvature____(r, s)
