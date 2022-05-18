# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 3:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.segments.segment.main import _2nCSCG_Segment
from root.config.main import rAnk, cOmm, sIze


class _2dCSCG_MeshLv0TraceElementsDo(FrozenOnly):
    """"""

    def __init__(self, lv0trace_elements):
        """"""
        self._elements_ = lv0trace_elements
        self._freeze_self_()


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

    def _Pr_update(self):
        """Update the information of segments on all local level-0-trace-elements.

        We do it in this way because we need to communication between cores. So we manually do it.
        """
        Segments = dict()
        for i in self._elements_:
            Segments[i] = list()

        mesh = self._elements_._mesh_
        lv0trace_elements = mesh.Lv0trace.elements
        MAP = mesh.cscg.trace.elements.map
        for i in mesh.base_cells:
            M = MAP[i]
            cell = mesh(i)

            _segments = dict()
            for Mi in M:
                _segments[Mi] = tuple()

            for j in cell:
                sub_cell = mesh(j)

                if sub_cell.IS.attached_to_Lv0cell_boundary:
                    if sub_cell.IS.attached_to_Lv0cell_U_boundary:
                        te = lv0trace_elements[M[0]]
                        xSignature = None
                        ySignature = self.___Pr_segment_finder___(sub_cell, 'U')
                        _segments[te.i] += ([xSignature, ySignature],)

                    if sub_cell.IS.attached_to_Lv0cell_D_boundary:
                        te = lv0trace_elements[M[1]]
                        xSignature = None
                        ySignature = self.___Pr_segment_finder___(sub_cell, 'D')
                        _segments[te.i] += ([xSignature, ySignature],)

                    if sub_cell.IS.attached_to_Lv0cell_L_boundary:
                        te = lv0trace_elements[M[2]]
                        xSignature = self.___Pr_segment_finder___(sub_cell, 'L')
                        ySignature = None
                        _segments[te.i] += ([xSignature, ySignature],)

                    if sub_cell.IS.attached_to_Lv0cell_R_boundary:
                        te = lv0trace_elements[M[3]]
                        xSignature = self.___Pr_segment_finder___(sub_cell, 'R')
                        ySignature = None
                        _segments[te.i] += ([xSignature, ySignature],)

                else:
                    pass

            for Mi in _segments:
                # noinspection PyTypeChecker,PyUnresolvedReferences
                Segments[Mi].append(_segments[Mi])

        ATAd = [dict() if _ != rAnk else None for _ in range(sIze)]

        for i in Segments: #go through all local trace-elements (level-0-trace-elements)
            if len(Segments[i]) == 1: # this trace-element data may should be communicated.
                te = lv0trace_elements[i]._te_
                if te.IS.shared_by_cores: # not on the mesh-boundary, then communicate it.
                    sc = te.shared_with_core
                    ATAd[sc][i] = Segments[i][0]

        ATAd = cOmm.alltoall(ATAd)

        for d_ri in ATAd:
            if d_ri is not None:
                for tei in d_ri:
                    assert tei in Segments
                    Segments[tei].append(d_ri[tei])

        for i in Segments:
            te = lv0trace_elements[i]
            LEN = len(Segments[i])
            if LEN == 1:
                assert te._te_.IS.on_mesh_boundary
                # noinspection PyTypeChecker
                Segments[i] = [_2nCSCG_Segment(te, *_) for _ in Segments[i][0]]
            elif LEN == 2:
                SGs0 = list()
                SGs1 = list()
                # noinspection PyTypeChecker
                for Ss in Segments[i][0]:
                    SGs0.append(_2nCSCG_Segment(te, *Ss))
                # noinspection PyTypeChecker
                for Ss in Segments[i][1]:
                    SGs1.append(_2nCSCG_Segment(te, *Ss))

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

        self._elements_._segments_ = Segments







if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/Lv0_trace/elements/do.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(1000, refinement_intensity=0.8)

    # print(mesh.base_cells.internal_segments._2bud_Sgs_)
    # from root.read.main import read
    # mesh = read('test_mesh.mi')

    # elements = mesh.Lv0trace.elements
    # elements.do._Pr_update()