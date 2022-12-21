# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 6/7/2022 3:59 PM
"""

from components.freeze.base import FrozenOnly


class iR_Gathering_Matrix_DoFind(FrozenOnly):
    """"""

    def __init__(self, igm):
        """"""
        self._GM_ = igm
        self._freeze_self_()

    def elements_and_local_indices_of_dof(self, i):
        """Find the local mesh elements and local indices of global dof #i.

        When we find nothing, we return None, Otherwise, we return two lists of mesh elements, and
        local indices respectively.

        """
        root_cells = list()
        local_indices = list()

        for e in self._GM_:  # go through all local root-cells
            gv = self._GM_[e]  # get the local gathering vector in each local root-cells
            if i in gv:
                root_cells.append(e)
                local_indices.append(gv.index(i))

        if len(root_cells) == 0:
            return None
        else:
            return root_cells, local_indices

    def elements_and_local_indices_of_dofs(self, dofs):
        """"""
        root_cells = dict()
        local_indices = dict()

        for i in dofs:
            root_cells[i] = list()
            local_indices[i] = list()

        for e in self._GM_:  # go through all local root-cells
            gv = self._GM_[e]  # get the local gathering vector in each local root-cells
            for i in dofs:
                if i in gv:
                    root_cells[i].append(e)
                    local_indices[i].append(gv.index(i))

        return root_cells, local_indices
