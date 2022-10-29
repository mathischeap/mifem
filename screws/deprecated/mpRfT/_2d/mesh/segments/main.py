# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/22/2022 10:04 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.segments.bcW import mpRfT2_Mesh_Segments_bcW
from objects.mpRfT._2d.mesh.segments.num import mpRfT2_Mesh_Segments_Num
from objects.mpRfT._2d.mesh.segments.visualize import mpRfT2_Mesh_Segments_Visualize
from objects.mpRfT._2d.mesh.segments.Wds import mpRfT2_Mesh_SegmentWiseDataStructure
from root.config.main import COMM, RANK, MASTER_RANK





class mpRfT2_Mesh_Segments(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._bcW_ = mpRfT2_Mesh_Segments_bcW(self)
        self._num_ = mpRfT2_Mesh_Segments_Num(self)
        self._visualize_ = None
        self._segments_ = dict()
        for seg in self:
            self._segments_[seg.__repr__()] = seg
        self._freeze_self_()

    def __iter__(self):
        """Go through all local segments.

        Returns
        -------

        """
        #---- first go through all local segments.
        for i in self._mesh_.basic_cells.trace_elements:
            segments = self._mesh_.basic_cells.trace_segments[i]
            for seg in segments:
                yield seg

        #--- then go through all internal segments.
        for i in self._mesh_.basic_cells.internal_segments:
            segments = self._mesh_.basic_cells.internal_segments[i]
            for seg in segments:
                yield seg

    def __getitem__(self, sg_rp):
        """Return a local segment according to its repr."""
        return self._segments_[sg_rp]

    @property
    def bcW(self):
        """Basic-cell-wise segments."""
        return self._bcW_

    @property
    def num(self):
        return self._num_

    @property
    def visualization(self):
        """Basic-cell-wise segments."""
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_Mesh_Segments_Visualize(self)
        return self._visualize_

    def ___Pr_find_N_for_all_segments___(self, SNM):
        """"""
        #---- parse the neighbours of segments first ----------------------------------------
        self.___Pr_parse_neighbors_for_all_segments___()
        #====================================================================================

        if SNM == 1:
            self.___Pr_Method_1___()
        else:
            raise NotImplementedError(f"SNM={SNM}.")

    def ___Pr_Method_1___(self):
        """The method No.1

            - If two neighbors match each other, we set the N of the segment to the higher one.
            - If Two neighbors do not match each other, we set the N of the segment to that of the
                smaller root-cell.

        Returns
        -------

        """
        for rp in self._mesh_.rcfc:
            cell = self._mesh_[rp]
            frame = cell.frame
            for edge in frame:
                EF = frame[edge]
                for segment in EF:
                    if segment.N is None:

                        neighbors = segment.neighbors
                        nN = segment.neighbors_N
                        nei0, nei1 = neighbors

                        if isinstance(nei0, str):
                            pass
                        else:
                            nei0 = nei0.__repr__()

                        if isinstance(nei1, str):
                            pass
                        else:
                            nei1 = nei1.__repr__()

                        if '-' not in nei0 and '-' not in nei1:
                            # two basic root-cells
                            segment._N_ = max(nN)

                        elif '-' in nei0 and '-' not in nei1:
                            segment._N_ = nN[0]

                        elif '-' not in nei0 and '-' in nei1:
                            segment._N_ = nN[1]

                        else: # two non-basic root-cells.

                            lv0 = len(nei0.split('-')[1])
                            lv1 = len(nei1.split('-')[1])

                            if lv0 == lv1:
                                segment._N_ = max(nN)
                            elif lv0 > lv1:
                                segment._N_ = nN[0]
                            else:
                                segment._N_ = nN[1]

                    else:
                        pass

    def ___Pr_parse_neighbors_for_all_segments___(self):
        """"""
        #------------ first we find all cell neighbors -------------------------------------
        for rp in self._mesh_.rcfc:
            cell = self._mesh_[rp]
            frame = cell.frame
            for edge in frame:
                EF = frame[edge]
                for segment in EF:
                    segment._neighbors_ .append (cell)
                    segment._neighbors_N_.append(cell.N)

        #--------- Then we try to find the boundary neighbors ------------------------------
        cscg = self._mesh_.cscg
        RTE = cscg.boundaries.range_of_trace_elements

        for bn in RTE:
            trace_elements = RTE[bn]
            for i in trace_elements:
                segments = self._mesh_.basic_cells.trace_segments[i]
                for seg in segments:
                    seg._neighbors_ .append (bn)
                    assert len(seg._neighbors_) == 2, f"must be the case!"
                    assert len(seg._neighbors_N_) == 1, f"must be the case!"

        neighborsTbc = dict()
        NNTbc = dict()

        for seg in self:
            self.num._local_ += 1 # local number of segments.
            if len(seg.neighbors) == 1:
                assert len(seg._neighbors_N_) == 1, f"must be the case!"
                rp = seg.__repr__()
                neighborsTbc[rp] = seg.neighbors[0].__repr__()
                NNTbc[rp] = seg._neighbors_N_[0]

        neighborsTbc = COMM.gather(neighborsTbc, root=MASTER_RANK)
        NNTbc = COMM.gather(NNTbc, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            NEI = dict()
            NNN = dict()
            for i, B_ in enumerate(neighborsTbc):
                N_ = NNTbc[i]
                for key in B_:
                    if key in NEI:
                        NEI[key].append(B_[key])
                        NNN[key].append(N_[key])
                    else:
                        NEI[key] = [B_[key],]
                        NNN[key] = [N_[key],]
        else:
            NNN = None
            NEI = None

        NNN = COMM.bcast(NNN, root=MASTER_RANK)
        NEI = COMM.bcast(NEI, root=MASTER_RANK)

        for seg in self:
            if len(seg._neighbors_) == 1:
                rp = seg.__repr__()
                if rp in NEI:
                    c0, c1 = NEI[rp]
                    N0, N1 = NNN[rp]
                    c = seg._neighbors_[0].__repr__()
                    if c0 == c:
                        seg._neighbors_ .append(c1)
                        seg._neighbors_N_.append(N1)
                    elif c1 == c:
                        seg._neighbors_ .append (c0)
                        seg._neighbors_N_.append(N0)
                    else:
                        raise Exception()

    @property
    def Wds(self):
        return mpRfT2_Mesh_SegmentWiseDataStructure(self)



if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/segments/main.py

    # from objects.mpRfT._2d.master import MeshGenerator
    # mesh = MeshGenerator('rectangle')([3,3], 2, show_info=True)

    from __init__ import rfT2
    mesh = rfT2.rm(10)

    # for seg in mesh.segments:
    #     print(seg.neighbors, seg._neighbors_N_)
    # segments = mesh.segments
    # bcW = mesh.segments.bcW
    #
    # for sg in segments:
    #     print(sg, sg.N)
