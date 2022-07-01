# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/13 9:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT.base.forms.base.partial.dofs.include import mpRfT_Form_PartialDofs_Include
from objects.mpRfT.base.forms.base.partial.dofs.compile import mpRfT_Form_PartialDofs_Compile


class PartialDofs(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        assert 'mpRfT_form' in f.standard_properties.tags
        self._f_ = f
        self._mesh_ = f.mesh
        self._indicators_ = dict()
        # keys are root-cells, values are the indicator of the local dofs if the root-cells.
        self._include_ = mpRfT_Form_PartialDofs_Include(self)
        self._compile_ = mpRfT_Form_PartialDofs_Compile(self)
        self._freeze_self_()

    @property
    def include(self):
        return self._include_

    @property
    def compile(self):
        return self._compile_






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/base/forms/base/partial/dofs/main.py
    pass
