# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/28/2022 12:24 PM
"""
from objects.CSCG.base.discrete_fields.base.main import CSCG_DiscreteField


class _3dCSCG_DiscreteField(CSCG_DiscreteField):
    """"""

    def __init__(self, mesh, coordinates, values, name, structured=False, grid=None):
        """"""
        super(_3dCSCG_DiscreteField, self).__init__(mesh, coordinates, values, name,
                                                    structured=structured, grid=grid)
