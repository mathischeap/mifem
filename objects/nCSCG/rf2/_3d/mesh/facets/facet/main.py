# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12:23 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _3nCSCG_Facet(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        self._freeze_self_()


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
