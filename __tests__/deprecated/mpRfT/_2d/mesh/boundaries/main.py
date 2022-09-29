# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/30 5:54 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_Mesh_Boundaries(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._names_ = mesh.cscg.boundaries.names
        self._rrc_ = None
        self._freeze_self_()

    @property
    def names(self):
        return self._names_

    @property
    def rrc(self):
        """range of root cells.

        The local root cells and their edges on each mesh boundary.
        """
        if self._rrc_ is None:
            ROE = self._mesh_.cscg.boundaries.range_of_element_edges
            rrc = dict()
            for bn in ROE:
                elements = ROE[bn]
                rrc_bn = list()
                for element_edge in elements:
                    element, edge = element_edge[:-1], element_edge[-1]
                    element = int(element)

                    ATTR = f"attached_to_basic_cell_{edge}_boundary"

                    basic_cell = self._mesh_[element]
                    for root_cells_index in basic_cell:
                        rc = self._mesh_[root_cells_index]

                        if rc.IS.attached_to_basic_cell_boundary:

                            if getattr(rc.IS, ATTR):
                                rrc_bn.append((rc.__repr__(), edge))
                            else:
                                pass

                        else:
                            pass

                rrc[bn] = rrc_bn

            self._rrc_ = rrc

        return self._rrc_






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/boundaries/main.py

    from __init__ import rfT2

    mesh = rfT2.rm(100, mesh_pool='crazy')

    print(mesh.boundaries.rrc)
