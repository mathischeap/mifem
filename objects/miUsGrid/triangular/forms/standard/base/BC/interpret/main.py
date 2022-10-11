# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/7/2022 7:25 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from objects.miUsGrid.triangular.forms.standard.base.BC.interpret.local import miUsTriangle_SF_BC_Interpret_Local

class miUsTriangle_SF_BC_Interpret(FrozenOnly):
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
            self._local_ = miUsTriangle_SF_BC_Interpret_Local(self._sf_)
        return self._local_






if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/base/BC/interpret/main.py
    pass
