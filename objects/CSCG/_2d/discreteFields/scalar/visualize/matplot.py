# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/30/2022 11:50 PM
"""
from components.freeze.base import FrozenOnly
from components.matplotWrappers.contour import contour, contourf


class _2cCSCG_DS_VisualizeMatplot(FrozenOnly):
    """"""

    def __init__(self, ds):
        """"""
        self._ds_ = ds
        self._freeze_self_()

    def __call__(self, **kwargs):
        return self.contour(**kwargs)

    def contour(self, **kwargs):
        """

        Parameters
        ----------
        kwargs

        Returns
        -------

        """
        xy = self._ds_.coordinates
        v = self._ds_.values
        return contour(xy, v, **kwargs)

    def contourf(self, **kwargs):
        """"""
        xy = self._ds_.coordinates
        v = self._ds_.values
        return contourf(xy, v, **kwargs)
