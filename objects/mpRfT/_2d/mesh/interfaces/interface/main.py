# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/20 5:48 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.interfaces.interface.IS import mpRfT2_Mesh_InterfaceIS
from root.config.main import rAnk, sIze


class mpRfT2_Mesh_Interface(FrozenOnly):
    """"""

    def __init__(self, mesh, rp, frame, direction, omb=False, symmetric=False):
        """

        Parameters
        ----------
        mesh
        rp : str
            The repr of this interface.
        frame : str
            The frame info
        direction : str
            'UD' or 'LR'
        omb : bool
            If this interface is on mesh boundaries?
        symmetric : bool
            If the basic-cells on two are symmetric.
        """
        self._mesh_ = mesh
        self._rp = rp
        self._direction_ = direction
        self._frame_ = frame
        self._omb_ = omb
        self._symmetric_ = symmetric

        if omb:
            self._swc_ = None
            sg_rp = self.__repr__().split('>')[1]
            frame = self.frame
            assert isinstance(frame, str)
            cell, edge = frame[:-2], frame[-1]

            segments = self._mesh_[cell].frame[edge]
            for seg in segments:
                assert seg.__repr__() == sg_rp
            self._segments_ = (sg_rp,)
            self._LEN_ = 1
        else:
            if isinstance(frame, list) and len(frame) == 2 :
                cell0 = frame[0][:-2]
                cell1 = frame[1][:-2]
                assert cell0 in self._mesh_.rcfc

                seg0 = self._mesh_[cell0].frame[frame[0][-1]]
                assert len(seg0) == 1
                for seg in seg0:
                    assert seg is self._mesh_.segments[self.__repr__().split('>')[1]]

                self._segments_=(seg.__repr__(),)

                if cell1 in self._mesh_.rcfc:
                    self._swc_ = None

                    seg1 = self._mesh_[cell1].frame[frame[1][-1]]
                    assert len(seg1) == 1
                    for seg in seg1:
                        assert seg is self._mesh_.segments[self.__repr__().split('>')[1]]

                else:
                    if '-' in cell1:
                        cell1 = int(cell1.split('-')[0])
                    else:
                        cell1 = int(cell1)

                    core = self._mesh_.cscg.do.find.slave_of_element(cell1)
                    assert core != rAnk
                    self._swc_ = core

                self._LEN_ = 1

            elif isinstance(frame, str):
                cell, edge = frame[:-2], frame[-1]

                segS = self.__repr__().split('>')[1].split('|')
                self._LEN_ = len(segS)

                for seg in segS:
                    assert seg in self._mesh_.segments._segments_
                self._segments_ = tuple(segS)

                if cell not in self._mesh_.rcfc:
                    if '-' in cell:
                        cell = int(cell.split('-')[0])
                    else:
                        cell = int(cell)

                    core = self._mesh_.cscg.do.find.slave_of_element(cell)
                    assert core != rAnk
                    self._swc_ = core

                else:
                    segments = self._mesh_[cell].frame[edge]

                    assert len(segments) > 1
                    self._swc_ = -1
                    for i, seg in enumerate(segments):

                        SEG = self._mesh_.segments[segS[i]]
                        assert seg is SEG

                        neighbors = seg.neighbors
                        n0, n1 = neighbors
                        if isinstance(n1, str):
                            assert n0 is self._mesh_[cell]
                            if '-' in n1:
                                n1 = int(n1.split('-')[0])
                            else:
                                n1 = int(n1)

                            core = self._mesh_.cscg.do.find.slave_of_element(n1)
                            assert core != rAnk
                            if self._swc_ == -1:
                                self._swc_ = core
                            else:
                                assert self._swc_ == core
                        else:
                            assert n1.__repr__() in self._mesh_.rcfc
                            self._swc_ = None
            else:
                raise Exception()

        if self._swc_ is not None:
            assert self._swc_ in range(0, sIze) and self._swc_ != rAnk

        self._IS_ = mpRfT2_Mesh_InterfaceIS(self)
        self._freeze_self_()

    def __repr__(self):
        return self._rp


    def __len__(self):
        """How many segments are involved in this interface."""
        return self._LEN_


    @property
    def frame(self):
        """This interface is also this frame. This frame may not in this core."""
        return self._frame_

    @property
    def direction(self):
        return self._direction_

    @property
    def IS(self):
        return self._IS_

    @property
    def shared_with_core(self):
        """return True if this interface presents in two cores else False."""
        return self._swc_

    @property
    def segments(self):
        """Return the repr of the segments of this interface."""
        return self._segments_



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
