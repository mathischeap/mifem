# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 1:19 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class FacetsDo(FrozenOnly):
    """"""

    def __init__(self, facets):
        """"""
        self._facets_ = facets
        self._freeze_self_()

    def _Pr_update(self):
        """"""




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
