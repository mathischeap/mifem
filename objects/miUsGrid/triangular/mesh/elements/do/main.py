# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/19 9:52 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np
from objects.miUsGrid.triangular.mesh.elements.do.find import miUsGrid_TriangularMesh_Elements_DO_FIND


class miUsGrid_TriangularMesh_Elements_DO(FrozenOnly):
    """"""

    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._find_ = miUsGrid_TriangularMesh_Elements_DO_FIND(elements)
        self._freeze_self_()

    @property
    def find(self):
        return self._find_

    def generate_grid_data(self):
        """We generate the grid data using the coordinate transformation instead of the cells and
        points. We do this for checking the coordinate transformation.

        """
        x = np.array([-1, 1])
        y = np.array([-1, 1,])
        x, y = np.meshgrid(x, y, indexing='ij')

        GD = dict()

        for i in self._elements_:
            element = self._elements_[i]
            X, Y = element.coordinate_transformation.mapping(x, y)

            ex = np.concatenate([X[0,:], X[1,:][::-1]])
            ey = np.concatenate([Y[0,:], Y[1,:][::-1]])
            center = element.coordinate_transformation.mapping(0, 0)
            singular_vertex = element.coordinate_transformation.mapping(-1, -1)

            U_edge = element.coordinate_transformation.mapping(np.array([-0.9,-0.9]), np.array([-0.9,0.9]))
            D_edge = element.coordinate_transformation.mapping(np.array([0.9,0.9]), np.array([-0.9,0.9]))
            L_edge = element.coordinate_transformation.mapping(np.array([-0.9,0.9]), np.array([-0.9,-0.9]))
            R_edge = element.coordinate_transformation.mapping(np.array([-0.9,0.9]), np.array([0.9,0.9]))

            GD[i] = [(ex, ey), center, singular_vertex, [U_edge, D_edge, L_edge, R_edge]]

        return GD



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
