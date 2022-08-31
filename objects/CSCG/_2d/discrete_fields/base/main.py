# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/28/2022 12:24 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.CSCG.base.discrete_fields.base.main import CSCG_DiscreteField


class _2dCSCG_DiscreteField(CSCG_DiscreteField):
    """"""

    def __init__(self, mesh, coordinates, values, name, structured=False, linspaces=None):
        """"""
        super(_2dCSCG_DiscreteField, self).__init__(mesh, coordinates, values, name,  structured=structured, linspaces=linspaces)




if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/discrete_fields/base/main.py
    pass
