# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/25 2:09 PM

A form in space [H^1]^3 actually.

"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class classname(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
