# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/27 5:22 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly






class mpRfT2_Segment_IS(FrozenOnly):
    """"""

    def __init__(self, segment):
        """"""
        self._segment_ = segment
        self._omb_ = None
        self._freeze_self_()

    @property
    def on_mesh_boundary(self):
        if self._omb_ is None:
            nei = self._segment_.neighbors
            self._omb_ = any([_ in self._segment_.where._mesh_.boundaries.names for _ in nei])
        return self._omb_






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/segments/segment/IS.py
    pass
