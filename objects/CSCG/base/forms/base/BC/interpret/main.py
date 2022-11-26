# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/13 9:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
from objects.CSCG.base.forms.base.BC.interpret.local import CSCG_FORM_BC_Interpret_Local

class CSCG_FORM_BC_Interpret(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._mesh_ = sf.mesh
        self._local_ = None
        self._freeze_self_()

    def RESET_cache(self):
        self._local_ = None

    @property
    def local(self):
        """We interpret the BC in local elements."""
        if self._local_ is None:
            self._local_ = CSCG_FORM_BC_Interpret_Local(self._sf_)
        return self._local_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
