# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class classname(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        self._freeze_self_()


if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/form/trace/base/main.py
    pass
