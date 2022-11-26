# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/15 9:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_S1F_Numbering_Local(FrozenOnly):
    """"""

    def __init__(self, numbering):
        """"""
        self._numbering_ = numbering
        self._cache_ = dict()
        self._freeze_self_()

    def __getitem__(self, rc_rp):
        """"""
        N = self._numbering_._f_.N[rc_rp]
        if N in self._cache_:
            pass
        else:
            f = self._numbering_._f_
            if f.orientation == 'inner':
                self._cache_[N] =  getattr(self._numbering_._f_.mesh.space[N].local_numbering,
                                           '_2dCSCG_1Form_Inner')
            else:
                self._cache_[N] =  getattr(self._numbering_._f_.mesh.space[N].local_numbering,
                                           '_2dCSCG_1Form_Outer')
        return self._cache_[N]


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
