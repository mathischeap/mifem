# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 12:24 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.base.fields.base import FiledBase


class miUsGrid_FiledBase(FiledBase):

    def __init__(self, mesh, valid_time, name):
        super(miUsGrid_FiledBase, self).__init__(mesh, valid_time)
        self.standard_properties.name = name


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
