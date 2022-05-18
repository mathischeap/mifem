# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 1:01 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_SegmentsBCW(FrozenOnly):
    """"""
    def __init__(self, segments):
        """"""
        self._mesh_ = segments._mesh_
        self._segments_ = segments
        self._DICT_ = dict()
        self._freeze_self_()

    def __iter__(self):
        for i in self._mesh_.base_cells:
            yield i

    def __getitem__(self, i):
        """Return a dict, keys are __repr__ of B-C-W segments, values are the segments.

        Parameters
        ----------
        i : int
            The number of a lv0-cell (base-cell; cscg mesh element).

        Returns
        -------

        """
        if i in self._DICT_:
            pass
        else:
            Di = dict()
            internal_segments = self._mesh_.base_cells.internal_segments[i]
            for its in internal_segments:
                _r = its.__repr__()
                Di[_r] = its

            Tmap = self._mesh_.Lv0trace.elements.map[i]
            for t in Tmap:
                t_segments = self._mesh_.Lv0trace.elements.segments[t]
                for ts in t_segments:
                    _r = ts.__repr__()
                    Di[_r] = ts

            self._DICT_[i] = Di

        return self._DICT_[i]





if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/segments/BCW.py
    # from objects.nCSCG.rfT2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    # mesh = rm2(100, refinement_intensity=0.5)
    from root.read.main import read
    mesh = read('test_mesh.mi')

    segments = mesh.segments

    # for i in segments.BCW:
    #     print(i, segments.BCW[i])

    mesh.visualize()
    mesh.Lv0trace.visualize()
    mesh.segments.visualize()