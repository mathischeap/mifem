# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:31 PM
"""
from objects.base.form.base import FormBase
from objects.miUsGrid.base.form.base.BC.main import miUsGrid_Form_BC

class miUsGrid_FormBase(FormBase):
    """"""

    def __init__(self, mesh, space, name):
        """"""
        super(miUsGrid_FormBase, self).__init__(mesh, space, name)
        self.standard_properties.___PRIVATE_add_tag___('miUsGrid_form')

        self._discretize_ = None
        self._reconstruct_ = None
        self._cochain_ = None
        self._num_ = None
        self._BC_ = miUsGrid_Form_BC(self)

    @property
    def ndim(self):
        return self.mesh.ndim

    @property
    def reconstruct(self):
        return self._reconstruct_

    @property
    def discretize(self):
        return self._discretize_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def num(self):
        return self._num_