# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 7:46 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.refinements.future.main import mpRfT2_Mesh_FutureRefinements


class mpRfT2_Mesh_Refinements(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._current_ = None
        self._future_ = mpRfT2_Mesh_FutureRefinements(mesh)
        self._freeze_self_()

    def ___Pr_parse_current_refinements___(self):
        """"""
        refinements = dict()

        for ind in self._mesh_:
            cell = self._mesh_[ind]

            if cell.level > 0: # all level > 0 cells are recorded.
                rp = cell.__repr__()
                assert cell.N is not None, f"cell:[{rp}] N is None."
                refinements[rp] = cell.N

            else:
                if cell.N != self._mesh_.dN: # p-refined basic-root-cell, record it.
                    key = ind[0]
                    assert key in self._mesh_.cscg.elements and \
                           key not in refinements and \
                           cell.N is not None, f"must be!"
                    refinements[key] = cell.N
                else:
                    pass

        return refinements

    @property
    def current(self):
        """"""
        if self._current_ is None:
            self._current_ = self.___Pr_parse_current_refinements___()
        return self._current_

    @property
    def future(self):
        """The refinements for the evolved mesh."""
        return self._future_







if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
