# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/25/2022 1:38 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


class miUsGrid_TriangularMesh_Boundaries(FrozenOnly):
    """"""

    def __init__(self, mesh, conditions):
        """"""
        self._mesh_ = mesh
        self.___Pr_parse_boundaries___(conditions)
        self._freeze_self_()

    def ___Pr_parse_boundaries___(self, conditions):
        """"""
        assert isinstance(conditions, dict), f"Use dict to define boundaries pls."
        names = tuple(conditions.keys())
        for name in names:
            assert name.isalnum(), f"boundary name = {name} is invalid."
        self._names_ = names

        mesh = self._mesh_
        emp = mesh.elements.map
        self._range_element_edge = dict()
        for name in names:
            self._range_element_edge[name] = list()

        for c in emp:
            MAP = emp[c]
            for i, m in enumerate(MAP):
                if m == '': # the `i`th edge of element #c is empty. Must be a boundary.
                    element = mesh.elements[c]

                    if i == 0: # the 0th edge
                        vertex0 = element.coordinates[0]
                        vertex1 = element.coordinates[1]
                    elif i == 1: # the 1st edge
                        vertex0 = element.coordinates[0]
                        vertex1 = element.coordinates[2]
                    elif i == 2: # the 2nd edge
                        vertex0 = element.coordinates[1]
                        vertex1 = element.coordinates[2]
                    else:
                        raise Exception()

                    find = 0
                    for bn in conditions:
                        condition = conditions[bn]

                        if condition(*vertex0) and condition(*vertex1):
                            emp[c][i] = bn
                            find += 1

                    if find == 0: raise Exception(f"edge[{i}] of element[{c}] finds no suitable boundary.")
                    assert find == 1, f"edge[{i}] of element[{c}] finds more than one suitable boundary."

                    self._range_element_edge[emp[c][i]].append((c,i))

                else: #an internal edge, m[0] must be numeric. and must have `+` or `-`
                    assert m[0].isnumeric()
                    if '-' in m:
                        assert '+' not in m
                    else:
                        assert '+' in m

    @property
    def names(self):
        return self._names_

    @property
    def range_of_element_edge(self):
        """return a dict.

        Keys are all boundary names, values are the local (element number, edge index) on the corresponding
        boundary.
        """
        return self._range_element_edge



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/boundaries/main.py
    from objects.miUsGrid.triangular.__test__.Random.test_mesh import mesh


    print(mesh.boundaries.range_of_element_edge)
