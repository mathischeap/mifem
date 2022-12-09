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

    def __init__(self, f):
        """"""
        self._f_ = f
        self._mesh_ = f.mesh
        if f.BC.CF is not None:
            self._ct_ = f.BC.CF._current_time_
        else:
            self._ct_ = None
        self._local_ = None
        self._freeze_self_()

    @property
    def local(self):
        """We interpret the BC in local elements."""

        if self._local_ is None:
            # when there is not local interpret, we make it.
            self._local_ = CSCG_FORM_BC_Interpret_Local(self._f_)

        else:
            # the local interpret is already there.
            if self._f_.BC.CF._current_time_ != self._ct_:
                # if current_time has been changed, we clean the cochains.
                self._local_._cochains_ = None
            else:
                pass

        return self._local_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
