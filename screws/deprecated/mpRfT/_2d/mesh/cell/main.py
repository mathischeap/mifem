# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/04 9:58 PM
"""

import sys
if './' not in sys.path: sys.path.append('./')

from objects.mpRfT.base.mesh.cell import mpRfT_Mesh_Cell_Base
from objects.mpRfT._2d.mesh.cell.coordinate_transformation import mpRfT2_Mesh_Cell_CT
from objects.mpRfT._2d.mesh.cell.do import mpRfT2_Mesh_Cell_Do
from objects.mpRfT._2d.mesh.cell.IS import mpRfT2_Mesh_Cell_IS
from objects.mpRfT._2d.mesh.cell.frame.main import mpFfT2_CellFrame




class mpRfT2_Mesh_Cell(mpRfT_Mesh_Cell_Base):
    """"""
    def __init__(self, mesh, level, indices):
        super(mpRfT2_Mesh_Cell, self).__init__(mesh, level, indices)
        self.___sub_cells___ = None # when initialized, it is a root cell.
        self._ct_ = None
        self._do_ = None
        self._IS_ = None
        self._type_wrt_metric_ = None
        self._N_ = None # CANNOT set to be dN of the mesh.
        self._space_ = None
        self._frame_ = None # CANNOT initialize it here because we will only do it for root-cells.
        self.___indices_metric_N_key___ = None
        self.___metric_N_key___ = None
        self.___N_key___ = None
        self._freeze_self_()

    def __getitem__(self, indices):
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
        if self.___isroot___: # I am a root cell.
            yield self.indices
        else: # I have sub-cells, go through all sub-root-cells.
            for i in range(4): # 8 for 3d
                sub_cell = self.sub_cells[i]
                for j in sub_cell:
                    yield j

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = mpRfT2_Mesh_Cell_CT(self)
        return self._ct_

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = mpRfT2_Mesh_Cell_Do(self)
        return self._do_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = mpRfT2_Mesh_Cell_IS(self)
        return self._IS_

    @property
    def type_wrt_metric(self):
        """type w.r.t. metric."""
        if self._type_wrt_metric_ is None:
            ELE_TYPE = self.mesh.cscg.elements[self.indices[0]].type_wrt_metric
            self._type_wrt_metric_ = \
                ELE_TYPE.___CLASSIFY_mpRfT2_CELL_of_origin_and_delta___(
                    self.coordinate_transformation.origin_and_delta)
        return self._type_wrt_metric_

    @property
    def __Pr_indices_metric_N_key___(self):
        """{str} : If chaotic, the str is numeric."""
        if self.___indices_metric_N_key___ is None:
            if isinstance(self.type_wrt_metric.mark, int):
                self.___indices_metric_N_key___ = str(self.type_wrt_metric.mark)
            else:
                if '-' in self.__repr__():
                    self.___indices_metric_N_key___ = self.type_wrt_metric.mark + str(self.N) + \
                                                      self.__repr__().split('-')[1]
                else:
                    self.___indices_metric_N_key___ = self.type_wrt_metric.mark + str(self.N)

        return self.___indices_metric_N_key___

    @property
    def ___Pr_metric_N_key___(self):
        """{str} : If chaotic, the str is numeric."""
        if self.___metric_N_key___ is None:
            if isinstance(self.type_wrt_metric.mark, int):
                self.___metric_N_key___ = str(self.type_wrt_metric.mark)
            else:
                self.___metric_N_key___ = self.type_wrt_metric.mark + str(self.N)
        return self.___metric_N_key___

    @property
    def ___Pr_N_key___(self):
        """{str} : If chaotic, the str is numeric."""
        if self.___N_key___ is None:
            self.___N_key___ = str(self.N)
        return self.___N_key___

    @property
    def N(self):
        """"""
        return self._N_

    @property
    def space(self):
        if self._space_ is None:
            self._space_ = self.mesh.space[self.N]
        return self._space_

    @property
    def frame(self):
        assert self.IS.root, f"Can only access the frame of a root-cell."
        if self._frame_ is None:
            self._frame_ = mpFfT2_CellFrame(self)
        return self._frame_





if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/cell/main.py
    pass