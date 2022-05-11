# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/09 12:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.segments.segment.main import _2nCSCG_Segment


class BaseCellsInternalSegmentsDo(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def _Pr_find_self_segments(self, indices):
        """"""
        o, d = self._mesh_.do.find.origin_and_delta(indices)
        ox, oy = o
        ex, ey = ox + d, oy + d
        where = self._mesh_(indices[0])
        cell = self._mesh_(indices)

        ySignature = ''
        for i in indices[1:]:
            if i in (0, 1):
                ySignature += '0'
            elif i in (2, 3):
                ySignature += '1'
            else:
                raise Exception()

        xSignature = ''
        for i in indices[1:]:
            if i in (0, 2):
                xSignature += '0'
            elif i in (1, 3):
                xSignature += '1'
            else:
                raise Exception()

        UDLR = list()

        if not cell.IS.attached_to_Lv0cell_U_boundary:
            UDLR.append(_2nCSCG_Segment(where, ox, ySignature))

        if not cell.IS.attached_to_Lv0cell_D_boundary:
            UDLR.append(_2nCSCG_Segment(where, ex, ySignature))

        if not cell.IS.attached_to_Lv0cell_L_boundary:
            UDLR.append(_2nCSCG_Segment(where, xSignature, oy))

        if not cell.IS.attached_to_Lv0cell_R_boundary:
            UDLR.append(_2nCSCG_Segment(where, xSignature, ey))

        return UDLR

    def _Pr_update(self):
        """"""
        mesh = self._mesh_
        base_cells = mesh.base_cells
        _2bud_Sgs_ = self._mesh_.base_cells.internal_segments._2bud_Sgs_
        for i in base_cells: # go through all local lv0-cells
            if _2bud_Sgs_[i]:
                lv0cell = base_cells[i]
                segments = list()
                if lv0cell.IS.root:
                    pass # no internal segments, so segments[i] = list()
                else:
                    sub_cells = dict() # keys are level + 1, values are indices
                    for ind in lv0cell:
                        LI = len(ind)
                        if LI not in sub_cells:
                            sub_cells[LI] = list()
                        sub_cells[LI].append(ind)

                    max_level = max(sub_cells.keys())

                    while max_level > 1:
                        if max_level in sub_cells:
                            indices = sub_cells[max_level]
                            for ind in indices:
                                UDLR = self._Pr_find_self_segments(ind)
                                for _ in UDLR:
                                    NEW = True
                                    for e in segments:
                                        if e in  _:
                                            NEW = False
                                            break
                                    if NEW:
                                        segments.append(_)
                        else:
                            pass
                        max_level -= 1


                self._mesh_.base_cells.internal_segments._segments_[i] = segments

                _2bud_Sgs_[i] = False







if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/base_cells/internal_segments/do.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(1000, refinement_intensity=0.5)

    # from root.read.main import read
    # mesh = read('test_mesh.mi')

    trace = mesh.Lv0trace
    trace.elements.do.update()
    # print(trace.elements.segments)
    # trace.visualize()
    # mesh.visualize()

    ISs = mesh.base_cells.internal_segments

    ISs.do.update()