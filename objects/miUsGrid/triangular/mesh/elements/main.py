# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/16 9:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.mesh.elements.do.main import miUsGrid_TriangularMesh_Elements_DO
from objects.miUsGrid.triangular.mesh.elements.num import miUsGrid_TriangularMesh_Elements_Num
from objects.miUsGrid.triangular.mesh.elements.element.main import miUsGrid_TriangularMesh_Element


class miUsGrid_TriangularMesh_Elements(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._range_ = None
        self._distributions_ = None
        self._map_ = None
        self._cells_ = None
        self._points_ = None
        self.__num_GLOBAL_points__ = None
        self._do_ = miUsGrid_TriangularMesh_Elements_DO(self)
        self._num_ = miUsGrid_TriangularMesh_Elements_Num(self)
        self._elements_ = dict()
        self._freeze_self_()

    @property
    def map(self):
        """The element map of all local cells."""
        return self._map_

    @property
    def distributions(self):
        """A list of shape `(sIze,)`.

        How many local cells in each core? distribution[i] means core #`i` have this many local
         cells.
         """
        return self._distributions_

    @property
    def range(self):
        """The local cells numbers are exactly elements of this range."""
        return self._range_

    @property
    def indices(self):
        """Same as range."""
        return self._range_

    @property
    def cells(self):
        """A dict whose keys are local cells, values are corresponding vertexes (points)."""
        return self._cells_

    @property
    def points(self):
        """A dict, keys are local points (all the vertexes), values are corresponding coordinates."""
        return self._points_

    def __iter__(self):
        """Go through all local cells."""
        for i in self.range:
            yield i

    def __contains__(self, item):
        """If item is a local cell number?"""
        return item in self.range

    def __len__(self):
        """How many local cells?"""
        return len(self.range)

    def __getitem__(self, cell_number):
        """"""
        if cell_number not in self._elements_:
            self._elements_[cell_number] = miUsGrid_TriangularMesh_Element(self, cell_number)
        return self._elements_[cell_number]

    @property
    def do(self):
        return self._do_

    @property
    def num(self):
        return self._num_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
