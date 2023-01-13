# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:30 PM
"""
from objects.miUsGrid.base.form.base.main import miUsGrid_FormBase


class miUsTriangular_FormBase(miUsGrid_FormBase):
    """"""

    def __init__(self, mesh, space, name):
        """"""
        super(miUsTriangular_FormBase, self).__init__(mesh, space, name)
        self.standard_properties.___PRIVATE_add_tag___('miUsGrid_triangular_form')
