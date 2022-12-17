# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/02 1:59 PM
"""
from objects.CSCG._3d.discreteDields.base.main import _3dCSCG_DiscreteField



class _3dCSCG_DF_Scalar(_3dCSCG_DiscreteField):
    """Region wise scalar data."""

    def __init__(self, mesh, coordinates, values, name='no-name', structured=False, grid=None):
        """"""
        super(_3dCSCG_DF_Scalar, self).__init__(mesh, coordinates, values, name,
                                                structured=structured, grid=grid)
        assert self.vdim == 1, f"vdim must be 1, now it is {self.vdim}."
