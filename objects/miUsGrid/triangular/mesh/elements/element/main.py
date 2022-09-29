# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/19 11:33 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.mesh.elements.element.coordinate_transformation.main import miUsGrid_TriangularMesh_Element_CT
from objects.miUsGrid.triangular.mesh.elements.element.visualize import miUsGrid_TriangularMesh_Element_Visualize


class miUsGrid_TriangularMesh_Element(FrozenOnly):
    """"""

    def __init__(self, elements, i):
        """

        Parameters
        ----------
        elements
        i : int
            I am the #i element.

        """
        self._elements_ = elements
        self._i_ = i
        self._ct_ = None
        self._visualize_ = None
        self._freeze_self_()

    @property
    def i(self):
        """I am the `ith` element (triangular cell)."""
        return self._i_

    @property
    def coordinate_transformation(self):
        """"""
        if self._ct_ is None:
            self._ct_ = miUsGrid_TriangularMesh_Element_CT(self)
        return self._ct_

    @property
    def vertexes(self):
        """The global label of the three vertexes as cells in a sequence of 0->1->2."""
        return self._elements_.cells[self.i]

    @property
    def coordinates(self):
        """The coordinates of the three vertexes."""
        return [self._elements_.points[i] for i in self.vertexes]

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = miUsGrid_TriangularMesh_Element_Visualize(self)
        return self._visualize_



if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/elements/element/main.py
    pass
