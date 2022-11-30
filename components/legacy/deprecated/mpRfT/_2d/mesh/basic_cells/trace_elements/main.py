# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/22/2022 10:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from root.config.main import RANK, SIZE, COMM
from objects.mpRfT._2d.mesh.basic_cells.trace_elements.element import mpRfT2_Mesh_BasicCells_TraceElement
from objects.mpRfT._2d.mesh.basic_cells.trace_elements.visualize import mpRfT2_Mesh_BasicCells_TraceElements_Visualize
from objects.mpRfT._2d.mesh.segments.segment.main import mpRfT2_Segment


class mpRfT2_Mesh_BasicCells_TraceElements(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        COMM.barrier() #  because it will need global communications.
                       # Just make sure all cores reach this place.
        self._mesh_ = mesh
        self._elements_ = dict()
        self._visualize_ = None
        self.___Pr_parse_trace_elements_segments___()
        self._freeze_self_()

    @property
    def map(self):
        """dict : keys are level-0 cells, values are level-0-trace-elements around them."""
        return self._mesh_.cscg.trace.elements.map

    def __iter__(self):
        """Go through all local level-0-trace-elements"""
        for i in self._mesh_.cscg.trace.elements:
            yield i

    def __getitem__(self, i):
        if i in self._elements_:
            pass
        else:
            self._elements_[i] = mpRfT2_Mesh_BasicCells_TraceElement(self, self._mesh_, i)

        return self._elements_[i]

    def __contains__(self, i):
        return i in self._mesh_.cscg.trace.elements

    @property
    def segments(self):
        return self._segments_

    @staticmethod
    def ___Pr_segment_finder___(sub_cell, side):
        """We find the segment of a root-cell from itself.

        So this segment may not be valid because it can be divided from the other side!.
        """
        indices = sub_cell.indices[1:]
        signature = ''
        if side in 'UD':
            for i in indices:
                if i in (0, 1):
                    signature += '0'
                elif i in (2, 3):
                    signature += '1'
                else:
                    raise Exception()
        elif side in 'LR':
            for i in indices:
                if i in (0, 2):
                    signature += '0'
                elif i in (1, 3):
                    signature += '1'
                else:
                    raise Exception()
        else:
            raise Exception()
        return signature

    def ___Pr_parse_trace_elements_segments___(self):
        """Update the information of segments on all local level-0-trace-elements.

        We do it in this way because we need to communication between cores. So we manually do it.
        """
        Segments = dict()
        for i in self:
            Segments[i] = list()

        mesh = self._mesh_
        MAP = mesh.cscg.trace.elements.map
        for i in mesh.basic_cells:
            M = MAP[i]
            cell = mesh[i]

            _segments = dict()
            for Mi in M:
                _segments[Mi] = tuple()

            for j in cell:
                sub_cell = mesh[j]

                if sub_cell.whether.attached_to_basic_cell_boundary:
                    if sub_cell.whether.attached_to_basic_cell_U_boundary:
                        te = self[M[0]]
                        xSignature = None
                        ySignature = self.___Pr_segment_finder___(sub_cell, 'U')
                        _segments[te.i] += ([xSignature, ySignature],)

                    if sub_cell.whether.attached_to_basic_cell_D_boundary:
                        te = self[M[1]]
                        xSignature = None
                        ySignature = self.___Pr_segment_finder___(sub_cell, 'D')
                        _segments[te.i] += ([xSignature, ySignature],)

                    if sub_cell.whether.attached_to_basic_cell_L_boundary:
                        te = self[M[2]]
                        xSignature = self.___Pr_segment_finder___(sub_cell, 'L')
                        ySignature = None
                        _segments[te.i] += ([xSignature, ySignature],)

                    if sub_cell.whether.attached_to_basic_cell_R_boundary:
                        te = self[M[3]]
                        xSignature = self.___Pr_segment_finder___(sub_cell, 'R')
                        ySignature = None
                        _segments[te.i] += ([xSignature, ySignature],)

                else:
                    pass

            for Mi in _segments:
                # noinspection PyTypeChecker,PyUnresolvedReferences
                Segments[Mi].append(_segments[Mi])



        ATAd = [dict() if _ != RANK else None for _ in range(SIZE)]

        for i in Segments: #go through all local trace-elements (level-0-trace-elements)
            if len(Segments[i]) == 1: # these trace-element data may should be communicated.
                te = self[i]._te_
                if te.whether.shared_by_cores: # not on the mesh-boundary, then communicate it.
                    sc = te.shared_with_core
                    ATAd[sc][i] = Segments[i][0]

        ATAd = COMM.alltoall(ATAd)

        for d_ri in ATAd:
            if d_ri is not None:
                for tei in d_ri:
                    assert tei in Segments
                    Segments[tei].append(d_ri[tei])

        for i in Segments:
            te = self[i]
            LEN = len(Segments[i])
            if LEN == 1:
                assert te._te_.whether.on_mesh_boundary
                # noinspection PyTypeChecker
                Segments[i] = [mpRfT2_Segment(te, *_) for _ in Segments[i][0]]
            elif LEN == 2:
                SGs0 = list()
                SGs1 = list()
                # noinspection PyTypeChecker
                for Ss in Segments[i][0]:
                    SGs0.append(mpRfT2_Segment(te, *Ss))
                # noinspection PyTypeChecker
                for Ss in Segments[i][1]:
                    SGs1.append(mpRfT2_Segment(te, *Ss))

                SGs = list()
                for sg in SGs0:
                    IN = False
                    for _ in SGs1:
                        if sg in _:
                            IN = True
                            break
                    if IN:
                        SGs.append(sg)

                for sg in SGs1:
                    IN = False
                    for _ in SGs0:
                        if sg == _:
                            pass
                        elif sg in  _:
                            IN = True
                            break
                    if IN:
                        SGs.append(sg)

                Segments[i] = SGs

            else:
                raise Exception()

        self._segments_ = Segments

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_Mesh_BasicCells_TraceElements_Visualize(self)
        return self._visualize_




if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/basic_cells/trace_elements/main.py

    # from objects.mpRfT._2d.master import MeshGenerator
    # mesh = MeshGenerator('rectangle')([3,3], 2, show_info=True)

    from __init__ import rfT2
    mesh = rfT2.rm(50)
    mesh.basic_cells.trace_elements.visualize()
