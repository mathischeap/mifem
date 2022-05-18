# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/12 3:00 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.cell.space.do import _2nCSCG_MRF2_CellSpaceDo


class _2nCSCG_MeshCellSpace(FrozenOnly):
    """"""
    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._N_ = cell.mesh.space.dN # start with the default N.
        self._do_ = None
        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"{self._cell_.indices}-{self.N}:{self.body.__repr__()}"

    @property
    def N(self):
        """"""
        return self._N_

    @N.setter
    def N(self, N):
        """"""
        if N in self._cell_.mesh.space:
            self._N_ = N
        else:
            raise Exception(f'N={N} is not a valid space degree. Add the space first.')

    @property
    def body(self):
        """"""
        return self._cell_.mesh.space[self.N]

    @property
    def do(self):
        """"""
        if self._do_ is None:
            self._do_ = _2nCSCG_MRF2_CellSpaceDo(self)
        return self._do_

    @property
    def nodes(self):
        return self.body.nodes

    @property
    def GoN(self):
        return self._cell_.mesh.space.GoN[self.N]

    @property
    def GoN_ravel(self):
        return self._cell_.mesh.space.GoN_ravel[self.N]









if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/cell/space/main.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)

    for i in mesh:
        print(mesh(i).space)