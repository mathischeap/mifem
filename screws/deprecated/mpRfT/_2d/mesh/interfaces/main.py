# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/20 1:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.interfaces.interface.main import mpRfT2_Mesh_Interface
from objects.mpRfT._2d.mesh.interfaces.visualization import mpRfT2_Mesh_Interfaces_Visualization
from root.config.main import COMM, SIZE


class mpRfT2_Mesh_Interfaces(FrozenOnly):
    """An interface is a combination of segments who have only one basic-cell on, at least, one of
    its two sides."""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self.___Pr_parse_all_interfaces___()
        self._visualization_ = mpRfT2_Mesh_Interfaces_Visualization(self)
        self._freeze_self_()

    def ___Pr_parse_all_interfaces___(self):
        """"""
        IFs = set()

        mesh = self._mesh_
        BNs = mesh.boundaries.names
        Tmap = mesh.cscg.trace.elements.map
        T_elements = mesh.cscg.trace.elements

        local_further_seg = list()
        multi_segments_cOmm = list()
        local_further_seg_cOmm = list()

        for i in range(SIZE):
            multi_segments_cOmm.append(list())
            local_further_seg_cOmm.append(list())

        for rc_rp in mesh.rcfc:
            cell = mesh[rc_rp]
            frame = cell.frame
            for ei, edge in enumerate(frame):
                segments = frame[edge]
                LEN = len(segments)
                if LEN > 1:
                    # must be an interface
                    if_rp = segments.__repr__()
                    IFs.add(if_rp)

                    for seg in segments:

                        neighbors = seg.neighbors
                        if not isinstance(neighbors[1], str):
                            # two neighbors are both local, no need to communicate
                            pass
                        else:
                            # this if_rp will be communicated.
                            i0 = cell.indices[0]
                            TE = Tmap[i0][ei]
                            shared_with_core = T_elements[TE].shared_with_core
                            multi_segments_cOmm[shared_with_core].append(if_rp)

                        break # break after checking the first segment.

                elif LEN == 1:
                    for seg in segments: # get the only segment and name it seg

                        neighbors = seg.neighbors
                        if isinstance(neighbors[1], str):

                            assert neighbors[0] is cell, f"must be!"

                            if neighbors[1] in BNs:
                                IFs.add(seg.__repr__())

                            else:
                                # do nothing as we have to check the other side: it is also the only segment or not.
                                local_further_seg.append(seg.__repr__())

                                i0 = cell.indices[0]
                                TE = Tmap[i0][ei]
                                shared_with_core = T_elements[TE].shared_with_core
                                local_further_seg_cOmm[shared_with_core].append(seg.__repr__())

                        else:
                            # the two neighbors are both local
                            if cell is neighbors[0]:
                                other = neighbors[1]
                            else:
                                assert cell is neighbors[1], f"must be."
                                other = neighbors[0]

                            if edge == 'U':
                                other_edge = 'D'
                            elif edge == 'D':
                                other_edge = 'U'
                            elif edge == 'L':
                                other_edge = 'R'
                            else:
                                other_edge = 'L'

                            other_segments = other.frame[other_edge]

                            if len(other_segments) > 1:
                                pass
                            else:

                                sg_rp = seg.__repr__()
                                if sg_rp in IFs:
                                    pass
                                else:
                                    IFs.add(sg_rp)

        multi_segments_cOmm = COMM.alltoall(multi_segments_cOmm)
        local_further_seg_cOmm = COMM.alltoall(local_further_seg_cOmm)

        ___ = list()
        for _ in multi_segments_cOmm:
            ___.extend(_)
        multi_segments_cOmm = ___

        ___ = list()
        for _ in local_further_seg_cOmm:
            ___.extend(_)
        local_further_seg_cOmm = ___

        for i, sg_rp in enumerate(local_further_seg):
            if sg_rp in local_further_seg_cOmm:
                assert local_further_seg_cOmm.count(sg_rp) == 1
                IFs.add(sg_rp)
            else:
                TF = False
                for frame in multi_segments_cOmm:
                    if sg_rp in frame:
                        TF = True
                        IFs.add(frame)
                        break
                assert TF

        rp_temp = 'IF>'
        self.___interfaces___ = dict()

        for IF in IFs:
            if 'Frame-' in IF:
                sg_rps = IF.split('Frame-')[1][2:]
                if_rp = rp_temp + sg_rps

                frame_info = IF.split('|')[0]
                frame_info = frame_info.split(':Frame-')
                frame_info = frame_info[0] + '>' + frame_info[1]

                if frame_info[-1] in ('U', 'D'):
                    direction = 'LR'
                else:
                    direction = 'UD'
                interface = mpRfT2_Mesh_Interface(
                    self._mesh_, if_rp, frame_info, direction, omb=False, symmetric=False)
                self.___interfaces___[if_rp] = interface

            else:
                if_rp = rp_temp + IF
                seg = mesh.segments[IF]

                direction = seg.direction

                neighbors = seg.neighbors
                if isinstance(neighbors[1], str):
                    if neighbors[1] in BNs:
                        omb = True
                    else:
                        omb = False
                else:
                    omb = False

                if omb:
                    symmetric = False
                else:
                    symmetric = True

                frame_info = list()
                for nei in neighbors:

                    if not isinstance(nei, str):
                        frame = nei.frame

                        COUNT = 0
                        if direction == 'UD':
                            for edge in 'LR':
                                segments = frame[edge]
                                if segments.__repr__().split('|')[1] == IF:

                                    _ = segments.__repr__().split('|')[0]
                                    _ = _.split(':Frame-')

                                    frame_info.append(_[0]+'>'+_[1])


                                    COUNT += 1
                                    break
                        else:
                            for edge in 'UD':
                                segments = frame[edge]
                                if segments.__repr__().split('|')[1] == IF:

                                    _ = segments.__repr__().split('|')[0]
                                    _ = _.split(':Frame-')

                                    frame_info.append(_[0]+'>'+_[1])
                                    COUNT += 1
                                    break
                        assert COUNT == 1

                    else:
                        if nei in BNs:
                            frame_info = frame_info[0]
                        else:

                            __ = frame_info[0][-1]
                            if __ == 'D':
                                _ = 'U'
                            elif __ == 'U':
                                _ = 'D'
                            elif __ == 'L':
                                _ = 'R'
                            elif __ == 'R':
                                _ = 'L'
                            else:
                                raise Exception()
                            frame_info.append(str(nei) + '>' + _)

                interface = mpRfT2_Mesh_Interface(
                    self._mesh_, if_rp, frame_info, direction, omb=omb, symmetric=symmetric)
                self.___interfaces___[if_rp] = interface


    def __getitem__(self, if_rp):
        """Return a local interface whose repr is if_rp."""
        return self.___interfaces___[if_rp]

    def __iter__(self):
        """Go through all local interfaces."""
        for if_rp in self.___interfaces___:
            yield if_rp

    def __len__(self):
        """How many local interfaces?"""
        return len(self.___interfaces___)

    def __contains__(self, if_rp):
        """Check if if_rp is the repr of a local interface?"""
        return if_rp in self.___interfaces___


    @property
    def visualization(self):
        return self._visualization_




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/interfaces/main.py

    from __init__ import rfT2

    # for i in range(10):
    mesh = rfT2.rm(100, refinement_intensity=0.5)
    #     IF = mesh.interfaces

    # from __init__ import rfT2
    # N = 3
    # K = 2
    # # rfd = None
    # rfd = {'0-0':N, '0-10':N, '0-11':N, '0-12':N, '0-13':N, '0-2':N, '0-30':N, '0-31':N, '0-32':N, '0-33':N,
    #        '1-1':N, '1-00':N, '1-01':N, '1-02':N, '1-03':N, '1-2':N, '1-3':N,}
    #
    # mesh = rfT2.mesh('crazy', c=0)([K, K], N, rfd=rfd)

    IF = mesh.interfaces

    IF.visualization()

    # from root.config.main import RANK
    # if RANK == 0:
    #     for rp in IF:
    #         print(IF[rp].segments)

    # for seg in mesh.segments:
    #     print(seg.__repr__())

    # mesh.visualization()