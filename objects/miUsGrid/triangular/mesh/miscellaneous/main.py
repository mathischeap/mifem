# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 4:02 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from root.config.main import COMM, RANK, MASTER_RANK
import numpy as np


class miUsGrid_TriangularMesh_Miscellaneous(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._quadratic_triangle_data_CACHE_ = None
        self._TQ_cache_ = None
        self._freeze_self_()

    def RESET_cache(self):
        self._quadratic_triangle_data_CACHE_ = None
        self._TQ_cache_ = None

    @property
    def quadratic_triangle_data(self):
        """Return connections (2d, `.ravel('C')` will be the connection), offsets, types of
         a vtk unstructured quadratic triangle grid.
        """
        if self._quadratic_triangle_data_CACHE_ is None: # we cache it as it is not very fast.

            element_map = self._mesh_.elements.map
            p = 2
            element_map = COMM.gather(element_map, root=MASTER_RANK)
            cells = COMM.gather(self._mesh_.elements.cells, root=MASTER_RANK)

            if RANK == MASTER_RANK:
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

                    else: # when p = 1, no upper edge internal dofs, just pass
                        pass

                    #- STEP 2: element internal dofs ----------------------------------------------
                    pass # as no cell internal dofs ----------------------------------------------------

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

                    else: # when p = 1, no down edge internal dofs, just pass
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

                    else: # when p = 1, no right edge internal dofs, just pass
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

                self._quadratic_triangle_data_CACHE_ = (
                    CONNECTIONS.ravel('C'), offsets, types, global_indices, gathering
                )

            else:
                global_indices = [1, 3, 5, 6, 7, 8]
                self._quadratic_triangle_data_CACHE_ = (
                    None, None, None, global_indices, None
                )

        return self._quadratic_triangle_data_CACHE_

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


    def triangle_and_quad_data(self, p):
        """

        Parameters
        ----------
        p

        Returns
        -------

        """
        if self._TQ_cache_ is None or self._TQ_cache_[0] != p: # we cache it as it is not very fast.

            assert p % 1 == 0 and p >= 1, f"p={p} invalid, must be a positive integer."
            element_map = self._mesh_.elements.map
            element_map = COMM.gather(element_map, root=MASTER_RANK)
            cells = COMM.gather(self._mesh_.elements.cells, root=MASTER_RANK)

            if RANK == MASTER_RANK:
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

                    else: # when p = 1, no upper edge internal dofs, just pass
                        pass

                    #- STEP 2: element internal dofs ----------------------------------------------
                    if p > 1:

                        numbering_c[1:-1, 1:-1] = np.arange(
                            current, current + (p - 1)**2).reshape([(p-1), (p-1)], order='F')
                        current += (p - 1) ** 2

                    else:
                        pass

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

                    else: # when p = 1, no down edge internal dofs, just pass
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

                    else: # when p = 1, no right edge internal dofs, just pass
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



                CONNECTIONS = np.zeros((len(numbering), p*3 + p*(p-1)*4), dtype=int)
                offsets = np.zeros(p**2*len(numbering), dtype=int)

                global_indices = [int(p/2),] + list(np.arange(p+1, (p+1)**2))
                gathering = np.zeros((len(numbering), len(global_indices)), dtype=int)
                ct = 0
                basic_vtk_cell_number = 0
                for e in numbering:

                    assert basic_vtk_cell_number == e * p**2

                    ne = numbering[e]

                    CONNECTION_e = np.zeros(p*3 + p*(p-1)*4, dtype=int)
                    for i in range(p): # the p triangles adjacent to the singular vertex.
                        CONNECTION_e[i*3 : (i+1)*3] = [ne[0,0], ne[i+1,1], ne[i,1]]

                        ct += 3
                        offsets[basic_vtk_cell_number] = ct
                        basic_vtk_cell_number += 1

                    if p > 1: # the p * (p-1) quad cells.
                        for k in range(p-1):
                            for i in range(p):
                                j = k + 1
                                m = i + k * p
                                ind = p * 3 + m * 4
                                CONNECTION_e[ind : ind+4] = [
                                    ne[i  , j  ],
                                    ne[i+1, j  ],
                                    ne[i+1, j+1],
                                    ne[i  , j+1]
                                ]

                                ct += 4
                                offsets[basic_vtk_cell_number] = ct
                                basic_vtk_cell_number += 1

                    CONNECTIONS[e, :] = CONNECTION_e
                    gathering[e, :] = ne.ravel('F')[global_indices]

                types = np.tile(np.array([5,]*p + [9,]*p*(p-1)), len(numbering))

                self._TQ_cache_ = (
                    p, CONNECTIONS.ravel('C'), offsets, types, global_indices, gathering
                )

            else:
                global_indices = [int(p/2),] + list(np.arange(p+1, (p+1)**2))
                self._TQ_cache_ = (
                    p, None, None, None, global_indices, None
                )

        return self._TQ_cache_[1:]


if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/miscellaneous/main.py
    from __init__ import miTri

    fc = miTri.form('rand0', 2)
    mesh = fc.mesh

    a = mesh.miscellaneous.triangle_and_quad_data(p=2)
