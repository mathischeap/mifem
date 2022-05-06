# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/04 9:58 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.mesh.cell.main import nCSCG_RF2_MeshCell
from objects.nCSCG.rf2._3d.mesh.cell.do import _3nCSCG_CellDo
from objects.nCSCG.rf2._3d.mesh.cell.IS import _3nCSCG_CellIS
from objects.nCSCG.rf2._3d.mesh.cell.coordinate_transformation import _3nCSCG_CellCT


class _3nCSCG_RF2_MeshCell(nCSCG_RF2_MeshCell):
    """"""
    def __init__(self, mesh, level, indices):
        super(_3nCSCG_RF2_MeshCell, self).__init__(mesh, level, indices)
        self.___sub_cells___ = None
        self._do_ = None
        self._IS_ = None
        self._CT_ = None
        self._type_wrt_metric_ = None
        self._freeze_self_()

    def __call__(self, indices):
        """"""
        if isinstance(indices, int):
            indices = (indices,)
        else:
            pass

        cells = self.sub_cells
        for _, i in enumerate(indices):
            if cells is None:
                raise Exception(f"cell[{self.indices + indices[:_]}] has no level: "
                                f">{_ + self.level + 1}< sub-cells.")
            else:
                cell = cells[i]

            cells = cell.sub_cells

        return cell

    def __iter__(self):
        """"""
        if self.sub_cells is None: # I am the smallest cell.
            yield self.indices
        else: # I have sub-cells, go through all sub-cells.
            for i in range(8):
                sub_cell = self.sub_cells[i]
                for j in sub_cell:
                    yield j

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _3nCSCG_CellDo(self)
        return self._do_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _3nCSCG_CellIS(self)
        return self._IS_

    @property
    def coordinate_transformation(self):
        if self._CT_ is None:
            self._CT_ = _3nCSCG_CellCT(self)
        return self._CT_

    @property
    def type_wrt_metric(self):
        """"""
        if self._type_wrt_metric_ is None:
            ELE_TYPE = self.mesh.cscg.elements[self.indices[0]].type_wrt_metric
            self._type_wrt_metric_ = \
                ELE_TYPE.___CLASSIFY_3nCSCG_RF2_CELL_of_origin_and_delta___(
                    self.coordinate_transformation.origin_and_delta)
        return self._type_wrt_metric_




if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_3d/mesh/cell/main.py
    from objects.nCSCG.rf2._3d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None)
    if 0 in mesh.cscg.elements:
        c0 = mesh(0)
        # print(c0.indices, c0.level, c0.___sub_cells___)

        c0.do.refine()
        # c00 = mesh[(0,0)]
        # print(c00.indices, c00.level, c00.___sub_cells___)
        #
        c00  = c0(0)

        # print(c00 is C00)
        c00.do.refine()
        c000 = mesh((0,0,7))
        C000 = c0((0,7))
        print(c000 is C000, c00.IS.refined, C000.IS.refined)

        print(C000.type_wrt_metric.mark)
