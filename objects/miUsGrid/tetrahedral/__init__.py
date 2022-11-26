# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/04 10:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly


class classname(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
