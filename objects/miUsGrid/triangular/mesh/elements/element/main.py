# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/19 11:33 AM
"""
import sys

import numpy as np

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.mesh.elements.element.coordinate_transformation.main import miUsGrid_TriangularMesh_Element_CT
from objects.miUsGrid.triangular.mesh.elements.element.visualize import miUsGrid_TriangularMesh_Element_Visualize
from screws.functions._2d_space.distance import distance
from screws.functions._2d_space.angles_of_triangle import angles_of_triangle

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
    def edge_lengths(self):
        """return three floats referring to the length of the three edges; edge0: 0->1 (D edge),
        edge1: 0->2 (U edge), edge2: 2->1 (R edge).
        """
        p0, p1, p2 = self.coordinates
        len0 = distance(p0, p1)
        len1 = distance(p0, p2)
        len2 = distance(p2, p1)
        return len0, len1, len2

    @property
    def area(self):
        """The area of this triangle element."""
        a, b, c = self.edge_lengths
        s = (a + b + c) / 2
        return (s * (s - a) * (s - b) * (s - c)) ** 0.5

    @property
    def angles(self):
        """return the angles in degrees (0-180)."""
        return angles_of_triangle(*self.coordinates)

    @property
    def quality(self):
        """1: best, 0, worst.

        A best triangle has three 60 degree angles.
        A worst one has a 180 angle plus two 0 degree ones.
        """

        diff = sum([abs(_ - 60) for _ in self.angles])
        return  abs(1 - diff / 240)

    @property
    def coordinates_2d(self):
        """Return a np.array of shape (3,2). `[0,:]` refers to (x,y) of vertex 0. `[:,0]` refers to
        x coordinates, `[:,1]` refers to y coordinates.
        """
        return np.vstack(self.coordinates)

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = miUsGrid_TriangularMesh_Element_Visualize(self)
        return self._visualize_

    @property
    def shape_indicator(self):
        """If two elements have the same shape_indicator, their shapes are the same and metric will
        be same as well.
        """
        c2d = self.coordinates_2d
        c2d[:,0] -= c2d[0,0]
        c2d[:,1] -= c2d[0,1]
        c2d = ''.join(['%.4f'%i for i in c2d[1:,:].ravel('F')])
        return c2d




if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/elements/element/main.py
    from __init__ import miTri
    fc = miTri.form('st2', 2)
    mesh = fc.mesh

    for e in mesh.elements:
        element = mesh.elements[e]
        print(element.area)
