# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 1:42 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_SegmentsDo(FrozenOnly):
    """"""

    def __init__(self, segments):
        """"""
        self._segments_ = segments
        self._freeze_self_()

    def _Pr_update(self):
        """update the segment information.

        This should be called only manually.
        """
        self._segments_._BCW_ = None


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
