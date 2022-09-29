# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/22/2022 10:08 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

class mpRfT2_Mesh_Segments_bcW(FrozenOnly):
    """Basic-Cell-Wise segments."""

    def __init__(self, segments):
        """"""
        self._mesh_ = segments._mesh_
        self._segments_ = segments
        self._DICT_ = dict()
        TMAP = self._mesh_.basic_cells.trace_elements.map

        for i in self._mesh_.basic_cells:

            Di = dict()
            internal_segments = self._mesh_.basic_cells.internal_segments[i]
            for its in internal_segments:
                _r = its.__repr__()
                Di[_r] = its

            Tmap = TMAP[i]

            for t in Tmap:
                t_segments = self._mesh_.basic_cells.trace_segments[t]
                for ts in t_segments:
                    _r = ts.__repr__()
                    Di[_r] = ts

            self._DICT_[i] = Di

        self._freeze_self_()

    def __iter__(self):
        for i in self._mesh_.basic_cells:
            yield i

    def __getitem__(self, i):
        """Return a dict, keys are __repr__ of basic-cell-wise segments, values are the segments in
        the corresponding basic-cells also stored in a dictionary.

        Parameters
        ----------
        i : int
            The number of the basic cell; cscg mesh element.

        Returns
        -------

        """
        return self._DICT_[i]

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
