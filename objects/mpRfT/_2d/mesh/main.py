# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/17 7:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT.base.mesh.main import mpRfT_MeshBase
from objects.mpRfT._2d.mesh.basic_cells.main import mpRfT2_Mesh_BasicCells
from objects.mpRfT._2d.mesh.do.main import mpRfT2_Mesh_Do
from objects.mpRfT._2d.mesh.visualize import mpRfT2_Mesh_Visualize

class mpRfT2_Mesh(mpRfT_MeshBase):
    """"""

    def __init__(self, cscg, dN, cells_dict):
        """"""
        super(mpRfT2_Mesh, self).__init__(cscg, dN, cells_dict)
        self._do_ = mpRfT2_Mesh_Do(self)
        self._visualize_ = None
        self._freeze_self_()

    def ___Pr_initialize_cells___(self, dN, cells_dict):
        """

        Parameters
        ----------
        dN
        cells_dict : dict
            {__repr__ : N, ...} : keys repr of a root-cell, N the space degree of the root-cell.

            Therefore, all sub-cells must be set `N`.

        Returns
        -------

        """
        self.___basic_cells___ = mpRfT2_Mesh_BasicCells(self)

        for i in self.basic_cells: self[i]._N_ = dN

        if cells_dict is None: cells_dict = list()

        for rp in cells_dict:

            if '-' not in rp: # customize a basic-root-cell.
                i = int(rp)
                assert i in self.basic_cells, f"cell={rp} not a valid basic cell."
                self.basic_cells[i]._N_ = cells_dict[rp]

            else: # customize a sub-root-cell. N must be given.
                i, indices = rp.split('-')
                indices = tuple([int(i),] + [int(_) for _ in indices])

                cells = self.basic_cells
                for _, i in enumerate(indices):
                    if cells is None:
                        cell.do.___Pr_h_refine___()
                        cells = cell.sub_cells
                    else:
                        pass
                    cell = cells[i]
                    cells = cell.sub_cells

                cell._N_ = cells_dict[rp]

    @property
    def do(self):
        return self._do_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_Mesh_Visualize(self)
        return self._visualize_





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/main.py
    # from objects.mpRfT._2d.master import MeshGenerator
    # mesh = MeshGenerator('rectangle')([3,3], 2, show_info=True)

    from __init__ import rfT2
    mesh = rfT2.rm(100)
    #
    mesh.visualize(color_space=True)