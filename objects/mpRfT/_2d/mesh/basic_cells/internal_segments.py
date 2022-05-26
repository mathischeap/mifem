# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/22/2022 11:03 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.segments.segment.main import mpRfT2_Segment

class mpRfT2_Mesh_BasicCells_InternalSegments(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._segments_ = dict()
        self.___Pr_parse_internal_segments___()
        self._freeze_self_()

    def __getitem__(self, item):
        """"""
        return self._segments_[item]

    def __iter__(self):
        """go through all local basic-cells (cscg mesh elements.)"""
        for i in self._mesh_.cscg.elements:
            yield i

    def __contains__(self, item):
        return item in self._mesh_.cscg.elements

    def ___Pr_find_self_segments___(self, indices):
        """"""
        o, d = self._mesh_.do.find.origin_and_delta(indices)
        ox, oy = o
        ex, ey = ox + d, oy + d
        where = self._mesh_[indices[0]]
        cell = self._mesh_[indices]

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
            UDLR.append(mpRfT2_Segment(where, ox, ySignature))

        if not cell.IS.attached_to_Lv0cell_D_boundary:
            UDLR.append(mpRfT2_Segment(where, ex, ySignature))

        if not cell.IS.attached_to_Lv0cell_L_boundary:
            UDLR.append(mpRfT2_Segment(where, xSignature, oy))

        if not cell.IS.attached_to_Lv0cell_R_boundary:
            UDLR.append(mpRfT2_Segment(where, xSignature, ey))

        return UDLR

    def ___Pr_parse_internal_segments___(self):
        """"""
        mesh = self._mesh_
        base_cells = mesh.basic_cells
        for i in base_cells: # go through all local lv0-cells
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
                            UDLR = self.___Pr_find_self_segments___(ind)
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

            self._segments_[i] = segments






if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
