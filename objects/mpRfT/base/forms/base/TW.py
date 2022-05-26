# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 3:32 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')
from screws.freeze.base import FrozenOnly


class mpRfT_FormTW(FrozenOnly):
    """All time-dependent properties of a form."""
    def __init__(self, f):
        """"""
        self._f_ = f
        self._current_time_ = None

        self._func_ = None
        self._freeze_self_()

    @property
    def current_time(self):
        ct = self.func.current_time
        #-- check if current_time are all same when we have multiple time-dependent properties.
        return ct

    @current_time.setter
    def current_time(self, current_time):
        assert isinstance(current_time, (int, float))
        #---- set all TW properties' current_time
        self.func.current_time = current_time

    @property
    def func(self):
        return self._func_

    @func.setter
    def func(self, func):
        """"""
        self._f_.___Pr_check_func___(func)
        self._func_ = func




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
