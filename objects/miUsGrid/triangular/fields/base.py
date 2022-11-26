# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 12:23 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.base.fields.base import miUsGrid_FiledBase

class miUsGrid_TriangularFieldBase(miUsGrid_FiledBase):
    """"""

    def __init__(self, mesh, valid_time, name):
        """"""
        super(miUsGrid_TriangularFieldBase, self).__init__(mesh, valid_time, name)
        self.standard_properties.___PRIVATE_add_tag___('miUsTri_field')



if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/fields/base.py
    pass
