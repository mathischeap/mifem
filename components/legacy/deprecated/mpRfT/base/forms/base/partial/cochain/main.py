# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/13 9:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from objects.mpRfT.base.forms.base.partial.dofs.main import PartialDofs
from objects.mpRfT.base.forms.base.partial.cochain.include import mpRfT_Form_PartialCochain_Include
from objects.mpRfT.base.forms.base.partial.cochain.compile import mpRfT_Form_PartialCochain_Compile


class PartialCochain(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        assert 'mpRfT_form' in f.standard_properties.tags
        self._f_ = f
        self._mesh_ = f.mesh
        self._pd_ = PartialDofs(f)
        self._lc_ = dict() # local cochain; keys are self._pd_.keys, values are the corresponding local-cochain.
        self._include_ = mpRfT_Form_PartialCochain_Include(self)
        self._compile_ = mpRfT_Form_PartialCochain_Compile(self)
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
