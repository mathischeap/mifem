# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/07 3:46 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class FrameSegments(FrozenOnly):
    """"""

    def __init__(self, cell, edge_name, segments):
        """

        Parameters
        ----------
        cell
        edge_name : str
            {'U', 'D', 'L', 'R'}
        """
        self._cell_ = cell
        self._edge_ = edge_name
        self._segments_ = segments
        self._freeze_self_()

    def __repr__(self):
        """"""
        return self._cell_.__repr__() + ':Frame-' + self._edge_

    @property
    def edge(self):
        return self._edge_

    def __iter__(self):
        """Go through all segments on this edge of a root-cell."""
        for _r in self._segments_:
            yield self._segments_[_r]







if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/cell/frame/segments/main.py
    pass
