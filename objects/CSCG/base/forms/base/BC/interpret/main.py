# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/13 9:19 PM
"""
from components.freeze.main import FrozenOnly
from objects.CSCG.base.forms.base.BC.interpret.local import CSCG_FORM_BC_Interpret_Local

class CSCG_FORM_BC_Interpret(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._mesh_ = f.mesh
        # this will make local.dofs during initializing. Do not make it jit.
        self._local_ = CSCG_FORM_BC_Interpret_Local(self._f_)
        self._freeze_self_()

    @property
    def local(self):
        """We interpret the BC in local elements."""
        return self._local_