# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/07 3:46 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_CellSegments(FrozenOnly):
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
        keys = list(segments.keys())
        keys.sort()
        d = dict()
        self._direction_ = None
        for key in keys:
            d[key] = segments[key]
            if self._direction_ is None:
                self._direction_ = segments[key].direction
            else:
                assert self._direction_ == segments[key].direction
        self._segments_ = d
        self._rp = None
        self._type_wrt_metric_ = None
        self._od_ = None
        self._freeze_self_()

    def __repr__(self):
        """"""
        if self._rp is None:
            rp = self._cell_.__repr__() + ':Frame-' + self._edge_
            for seg in self._segments_:
                rp += '|' + seg
            self._rp = rp
        return self._rp

    @property
    def edge(self):
        return self._edge_

    def __iter__(self):
        """Go through all segments on this edge of a root-cell."""
        for _r in self._segments_:
            yield self._segments_[_r]

    def __len__(self):
        return len(self._segments_)

    @property
    def type_wrt_metric(self):
        """type w.r.t. metric."""
        if self._type_wrt_metric_ is None:
            self._type_wrt_metric_ = ''
            for seg in self:
                mark = seg.type_wrt_metric.mark
                if isinstance(mark, int):
                    self._type_wrt_metric_ = id(self)
                    break
                else:
                    self._type_wrt_metric_ += mark
        return self._type_wrt_metric_


    @property
    def origin_and_delta(self):
        if self._od_ is None:
            O, D = None, list()
            RANGE = 0
            for seg in self:
                o, d = seg.origin_and_delta

                if O is None:
                    O = o

                    if isinstance(O, tuple):
                        if self.direction == 'UD':
                            Oo = O[0]
                        else:
                            Oo = O[1]
                    else:
                        Oo = O

                else:
                    if isinstance(o, tuple):
                        if self.direction == 'UD':
                            o = o[0]
                        else:
                            o = o[1]
                    else:
                        pass

                    assert sum(D)  == o - Oo

                D.append(d)

            D = sum(D)
            self._od_ = O, D

        return self._od_

    @property
    def direction(self):
        return self._direction_



if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/cell/frame/segments/main.py
    pass
