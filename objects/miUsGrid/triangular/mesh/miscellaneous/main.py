# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 4:02 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from root.config.main import cOmm, rAnk, mAster_rank
import numpy as np


class miUsGrid_TriangularMesh_Miscellaneous(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def quadratic_triangle_data(self):
        """Return connections (2d, `.ravel('C')` will be the connection), offsets, types of
         a vtk unstructured quadratic triangle grid.
        """
        element_map = self._mesh_.elements.map
        p = 2
        element_map = cOmm.gather(element_map, root=mAster_rank)
        cells = cOmm.gather(self._mesh_.elements.cells, root=mAster_rank)

        if rAnk == mAster_rank:
            ___ = dict()
            for em in element_map:
                ___.update(em)
            element_map = ___
            ___ = dict()
            for cs in cells:
                ___.update(cs)
            cells = ___

            current = 0
            numbering = dict()
            point_numbering = - np.ones(self._mesh_.elements.num.GLOBAL_points)

            for c in range(self._mesh_.elements.num.GLOBAL_cells):
                # we number all element by element according to the element number
                numbering_c = - np.ones((p+1, p+1), dtype=int)
                emp = element_map[c]
                vertexes = cells[c]

                #- STEP 0: numbering the singular vertex --------------------------------------
                singular_vertex = vertexes[0]
                if point_numbering[singular_vertex] == -1: # this point not numbered
                    numbering_c[:, 0] = current
                    point_numbering[singular_vertex] = current
                    current += 1
                else:
                    numbering_c[:, 0] = point_numbering[singular_vertex]

                #- STEP 1: numbering Upper internal (no-vertex) dofs -------------------------
                if p > 1:
                    what = emp[1]
                    if what.isalnum(): # boundary
                        assert what in self._mesh_.boundaries.names, f"trivial check"
                        numbering_c[0,1:-1] = np.arange(current, current+p-1)
                        current += p-1
                    else:
                        edge = what[-1]
                        sign = what[-2]
                        other = int(what[:-2])
                        if other > c: # neighbour is not numbered.
                            numbering_c[0, 1:-1] = np.arange(current, current + p - 1)
                            current += p-1
                        elif other < c:
                            numbering_c[0, 1:-1] = \
                                self.___Pr_find_0dofs_of_element_edge___(numbering[other], edge, sign)
                        else:
                            raise Exception(f"element[{c}] is neighbour to self.")

                else: # when p = 1, no internal dofs, just pass
                    pass

                #- STEP 2: element internal dofs ----------------------------------------------
                pass # as no internal dofs ----------------------------------------------------

                #- STEP 3 numbering Down internal (no-vertex) dofs ----------------------------
                if p > 1:
                    what = emp[0]
                    if what.isalnum(): # boundary
                        assert what in self._mesh_.boundaries.names, f"trivial check"
                        numbering_c[-1,1:-1] = np.arange(current, current+p-1)
                        current += p-1
                    else:
                        edge = what[-1]
                        sign = what[-2]
                        other = int(what[:-2])
                        if other > c: # neighbour is not numbered.
                            numbering_c[-1, 1:-1] = np.arange(current, current + p - 1)
                            current += p-1
                        elif other < c:
                            numbering_c[-1, 1:-1] = \
                                self.___Pr_find_0dofs_of_element_edge___(numbering[other], edge, sign)
                        else:
                            raise Exception(f"element[{c}] is neighbour to self.")

                else: # when p = 1, no internal dofs, just pass
                    pass

                #- STEP 4: numbering vertex 2 -------------------------------------------------
                vertex2 = vertexes[2]
                if point_numbering[vertex2] == -1: # this point not numbered
                    numbering_c[0, -1] = current
                    point_numbering[vertex2] = current
                    current += 1
                else:
                    numbering_c[0, -1] = point_numbering[vertex2]

                #- STEP 5 numbering Right internal (no-vertex) dofs ----------------------------
                if p > 1:
                    what = emp[2]
                    if what.isalnum(): # boundary
                        assert what in self._mesh_.boundaries.names, f"trivial check"
                        numbering_c[1:-1,-1] = np.arange(current, current+p-1)
                        current += p-1
                    else:
                        edge = what[-1]
                        sign = what[-2]
                        other = int(what[:-2])
                        if other > c: # neighbour is not numbered.
                            numbering_c[1:-1,-1] = np.arange(current, current + p - 1)
                            current += p-1
                        elif other < c:
                            numbering_c[1:-1,-1] = \
                                self.___Pr_find_0dofs_of_element_edge___(numbering[other], edge, sign)
                        else:
                            raise Exception(f"element[{c}] is neighbour to self.")

                else: # when p = 1, no internal dofs, just pass
                    pass

                #- STEP 6: numbering vertex 1 ---------------------------------------------------
                vertex1 = vertexes[1]
                if point_numbering[vertex1] == -1: # this point not numbered
                    numbering_c[-1, -1] = current
                    point_numbering[vertex1] = current
                    current += 1
                else:
                    numbering_c[-1, -1] = point_numbering[vertex1]

                #=================================================================================

                numbering[c] = numbering_c

            #_____________________________________________________________________________________
            CONNECTIONS = np.zeros((len(numbering), 6), dtype=int)

            local_indices = [0, 8, 6, 5, 7, 3]
            global_indices = [1, 3, 5, 6, 7, 8]
            gathering = np.zeros((len(numbering), 6), dtype=int)
            for e in numbering:
                Ne = numbering[e].ravel('F')
                CONNECTIONS[e, :] = Ne[local_indices]
                gathering[e, :] = Ne[global_indices]

            offsets = np.linspace(6, 6*len(numbering), len(numbering), dtype=int)

            types = 22

            return CONNECTIONS.ravel('C'), offsets, types, \
                   global_indices, gathering

        else:
            global_indices = [1, 3, 5, 6, 7, 8]
            return None, None, None, global_indices, None

    @staticmethod
    def ___Pr_find_0dofs_of_element_edge___(numbering, edge, sign):
        """"""
        if edge == '0': # Down edge
            dofs = numbering[-1,1:-1]
        elif edge == '1': # Upper edge
            dofs = numbering[0, 1:-1]
        elif edge == '2': # Right edge
            dofs = numbering[1:-1, -1]
        else:
            raise Exception()

        if sign == '+':
            return dofs
        elif sign == '-':
            return dofs[::-1]
        else:
            raise Exception()



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/vtk/main.py
    from __init__ import miTri

    fc = miTri.form('rand0', 2)
    mesh = fc.mesh

    a = mesh.miscellaneous.quadratic_triangle_data
