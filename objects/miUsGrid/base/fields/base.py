# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 12:24 PM
"""
from objects.base.fields.base import FiledBase


class miUsGrid_FiledBase(FiledBase):

    def __init__(self, mesh, valid_time, name):
        super(miUsGrid_FiledBase, self).__init__(mesh, valid_time)
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('miUs_scalar_field')
