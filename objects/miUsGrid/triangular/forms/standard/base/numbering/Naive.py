# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/25/2022 12:23 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from root.config.main import RANK, MASTER_RANK, COMM, np
from tools.linearAlgebra.gathering.vector import Gathering_Vector
from tools.linearAlgebra.gathering.regular.matrix.main import Gathering_Matrix

class Naive(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

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

    @property
    def miUsTriangular_S0F_Inner(self):
        return self.miUsTriangular_S0F_Outer

    @property
    def miUsTriangular_S0F_Outer(self):
        """"""
        element_map = self._sf_.mesh.elements.map
        p = self._sf_.space.p
        element_map = COMM.gather(element_map, root=MASTER_RANK)
        cells = COMM.gather(self._sf_.mesh.elements.cells, root=MASTER_RANK)

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
            point_numbering = - np.ones(self._sf_.mesh.elements.num.GLOBAL_points)

            for c in range(self._sf_.mesh.elements.num.GLOBAL_cells):
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
                        assert what in self._sf_.mesh.boundaries.names, f"trivial check"
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

                #- STEP 2: element internal dofs ---------------------------------------------
                if p > 1:

                    numbering_c[1:-1, 1:-1] = np.arange(current, current + (p - 1)**2).reshape([(p-1), (p-1)], order='F')
                    current += (p - 1) ** 2

                else:
                    pass

                #- STEP 3 numbering Down internal (no-vertex) dofs ----------------------------
                if p > 1:
                    what = emp[0]
                    if what.isalnum(): # boundary
                        assert what in self._sf_.mesh.boundaries.names, f"trivial check"
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

                #- STEP 4: numbering vertex 2 --------------------------------------
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
                        assert what in self._sf_.mesh.boundaries.names, f"trivial check"
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

                #- STEP 6: numbering vertex 1 --------------------------------------
                vertex1 = vertexes[1]
                if point_numbering[vertex1] == -1: # this point not numbered
                    numbering_c[-1, -1] = current
                    point_numbering[vertex1] = current
                    current += 1
                else:
                    numbering_c[-1, -1] = point_numbering[vertex1]

                #===================================================================

                numbering[c] = numbering_c

            numbering_distribution = list()
            distributions = self._sf_.mesh.elements.distributions
            for core, dis in enumerate(distributions):
                eMap = dict()
                for c in dis:
                    eMap[c] = numbering[c]
                numbering_distribution.append(eMap)

        else:
            numbering_distribution = None

        numbering_distribution = COMM.scatter(numbering_distribution, root=MASTER_RANK)
        GV_dict = dict()
        for c in numbering_distribution:
            full_vec = np.concatenate(([numbering_distribution[c][0,0],], numbering_distribution[c][:,1:].ravel('F')))
            GV_dict[c] = Gathering_Vector(c, full_vec)

        GM = Gathering_Matrix(GV_dict, mesh_type='miUsGrid_TriangularMesh')

        return GM

    @property
    def miUsTriangular_S2F_Inner(self):
        return self.miUsTriangular_S2F_Outer

    @property
    def miUsTriangular_S2F_Outer(self):
        """"""
        p = self._sf_.space.p
        GV_dict = dict()
        for c in self._sf_.mesh.elements.indices:
            GV_dict[c] = Gathering_Vector(c, range(c*p**2, (c+1)*p**2))
        GM = Gathering_Matrix(GV_dict, mesh_type='miUsGrid_TriangularMesh')
        return GM


    @staticmethod
    def ___Pr_find_1dofs_of_element_edge___(numbering, edge, sign):
        """"""
        if edge == '0': # Down edge
            dofs = numbering[-1,:]
        elif edge == '1': # Upper edge
            dofs = numbering[0, :]
        elif edge == '2': # Right edge
            dofs = numbering[:, -1]
        else:
            raise Exception()

        if sign == '+':
            return dofs
        elif sign == '-':
            return dofs[::-1]
        else:
            raise Exception()
    
    @property
    def miUsTriangular_S1F_Outer(self):

        element_map = self._sf_.mesh.elements.map
        p = self._sf_.space.p
        element_map = COMM.gather(element_map, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ___ = dict()
            for em in element_map:
                ___.update(em)
            element_map = ___

            current = 0
            numbering_DY = dict()
            numbering_DX = dict()

            for c in range(self._sf_.mesh.elements.num.GLOBAL_cells):
                # we number all element by element according to the element number
                numbering_dy = - np.ones((p+1, p), dtype=int)
                numbering_dx = - np.ones((p, p), dtype=int)
                emp = element_map[c]

                # STEP 0 --------- dy on Upper edge ----------------------
                what = emp[1]
                if what.isalnum(): # boundary
                    assert what in self._sf_.mesh.boundaries.names, f"trivial check"
                    numbering_dy[0,:] = np.arange(current, current+p)
                    current += p
                else:
                    edge = what[-1]
                    sign = what[-2]
                    other = int(what[:-2])
                    if other > c: # neighbour is not numbered.
                        numbering_dy[0, :] = np.arange(current, current + p)
                        current += p
                    elif other < c:

                        if edge in '01':
                            numbering_dy[0, :] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DY[other], edge, sign)
                        else:
                            numbering_dy[0, :] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DX[other], edge, sign)

                    else:
                        raise Exception(f"element[{c}] is neighbour to self.")

                #-- STEP 1 ------ internal dy edges -------------------------------------

                if p > 1:

                    numbering_dy[1:-1, :] = np.arange(current, current + p*(p-1)) .reshape([(p-1), p], order='F')
                    current += p*(p-1)

                else:
                    pass

                #-- STEP 2 --- dy on Down edge -----------------------------------------------
                what = emp[0]
                if what.isalnum(): # boundary
                    assert what in self._sf_.mesh.boundaries.names, f"trivial check"
                    numbering_dy[-1,:] = np.arange(current, current+p)
                    current += p
                else:
                    edge = what[-1]
                    sign = what[-2]
                    other = int(what[:-2])
                    if other > c: # neighbour is not numbered.
                        numbering_dy[-1, :] = np.arange(current, current + p)
                        current += p
                    elif other < c:

                        if edge in '01':
                            numbering_dy[-1, :] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DY[other], edge, sign)
                        else:
                            numbering_dy[-1, :] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DX[other], edge, sign)
                    else:
                        raise Exception(f"element[{c}] is neighbour to self.")

                #-- STEP 3 ------ internal dx edges -------------------------------------
                if p > 1:
                    numbering_dx[:, :-1] = np.arange(current, current + p*(p-1)).reshape([p, (p-1)], order='F')
                    current += p*(p-1)
                else:
                    pass

                #-- STEP 4 --- dx on Right edge -----------------------------------------------
                what = emp[2]
                if what.isalnum(): # boundary
                    assert what in self._sf_.mesh.boundaries.names, f"trivial check"
                    numbering_dx[:,-1] = np.arange(current, current+p)
                    current += p
                else:
                    edge = what[-1]
                    sign = what[-2]
                    other = int(what[:-2])

                    if other > c: # neighbour is not numbered.
                        numbering_dx[:, -1] = np.arange(current, current + p)
                        current += p

                    elif other < c:

                        if edge in '01':
                            numbering_dx[:, -1] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DY[other], edge, sign)
                        else:
                            numbering_dx[:, -1] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DX[other], edge, sign)

                    else:
                        raise Exception(f"element[{c}] is neighbour to self.")

                #===================================================================================
                numbering_DY[c] = numbering_dy
                numbering_DX[c] = numbering_dx

            distributions = self._sf_.mesh.elements.distributions

            DisDy = list()
            DisDx = list()

            for core, dis in enumerate(distributions):
                dis_dy = dict()
                dis_dx = dict()
                for c in dis:
                    dis_dy[c] = numbering_DY[c]
                    dis_dx[c] = numbering_DX[c]
                DisDy.append(dis_dy)
                DisDx.append(dis_dx)
        else:
            DisDy = None
            DisDx = None

        DisDy = COMM.scatter(DisDy, root=MASTER_RANK)
        DisDx = COMM.scatter(DisDx, root=MASTER_RANK)

        GV_dict = dict()
        for c in self._sf_.mesh.elements.indices:
            full_vec = np.concatenate([DisDy[c].ravel('F'), DisDx[c].ravel('F')])
            GV_dict[c] = Gathering_Vector(c, full_vec)

        GM = Gathering_Matrix(GV_dict, mesh_type='miUsGrid_TriangularMesh')

        return GM

    @property
    def miUsTriangular_S1F_Inner(self):

        element_map = self._sf_.mesh.elements.map
        p = self._sf_.space.p
        element_map = COMM.gather(element_map, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ___ = dict()
            for em in element_map:
                ___.update(em)
            element_map = ___

            current = 0
            numbering_DX = dict()
            numbering_DY = dict()

            for c in range(self._sf_.mesh.elements.num.GLOBAL_cells):
                # we number all element by element according to the element number
                numbering_dx = - np.ones((p, p), dtype=int)
                numbering_dy = - np.ones((p+1, p), dtype=int)
                emp = element_map[c]

                #-- STEP 0 ------ internal dx edges -------------------------------------
                if p > 1:
                    numbering_dx[:, :-1] = np.arange(current, current + p*(p-1)).reshape([p, (p-1)], order='F')
                    current += p*(p-1)
                else:
                    pass

                #-- STEP 1 --- dx on Right edge -----------------------------------------------
                what = emp[2]
                if what.isalnum(): # boundary
                    assert what in self._sf_.mesh.boundaries.names, f"trivial check"
                    numbering_dx[:,-1] = np.arange(current, current+p)
                    current += p
                else:
                    edge = what[-1]
                    sign = what[-2]
                    other = int(what[:-2])

                    if other > c: # neighbour is not numbered.
                        numbering_dx[:, -1] = np.arange(current, current + p)
                        current += p

                    elif other < c:

                        if edge in '01':
                            numbering_dx[:, -1] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DY[other], edge, sign)
                        else:
                            numbering_dx[:, -1] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DX[other], edge, sign)

                    else:
                        raise Exception(f"element[{c}] is neighbour to self.")
                # STEP 2 --------- dy on Upper edge ----------------------
                what = emp[1]
                if what.isalnum(): # boundary
                    assert what in self._sf_.mesh.boundaries.names, f"trivial check"
                    numbering_dy[0,:] = np.arange(current, current+p)
                    current += p
                else:
                    edge = what[-1]
                    sign = what[-2]
                    other = int(what[:-2])
                    if other > c: # neighbour is not numbered.
                        numbering_dy[0, :] = np.arange(current, current + p)
                        current += p
                    elif other < c:

                        if edge in '01':
                            numbering_dy[0, :] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DY[other], edge, sign)
                        else:
                            numbering_dy[0, :] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DX[other], edge, sign)

                    else:
                        raise Exception(f"element[{c}] is neighbour to self.")

                #-- STEP 3 ------ internal dy edges -------------------------------------
                if p > 1:
                    numbering_dy[1:-1, :] = np.arange(current, current + p*(p-1)) .reshape([(p-1), p], order='F')
                    current += p*(p-1)

                else:
                    pass

                #-- STEP 4 --- dy on Down edge -----------------------------------------------
                what = emp[0]
                if what.isalnum(): # boundary
                    assert what in self._sf_.mesh.boundaries.names, f"trivial check"
                    numbering_dy[-1,:] = np.arange(current, current+p)
                    current += p
                else:
                    edge = what[-1]
                    sign = what[-2]
                    other = int(what[:-2])
                    if other > c: # neighbour is not numbered.
                        numbering_dy[-1, :] = np.arange(current, current + p)
                        current += p
                    elif other < c:

                        if edge in '01':
                            numbering_dy[-1, :] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DY[other], edge, sign)
                        else:
                            numbering_dy[-1, :] = \
                                self.___Pr_find_1dofs_of_element_edge___(numbering_DX[other], edge, sign)
                    else:
                        raise Exception(f"element[{c}] is neighbour to self.")

                #===================================================================================
                numbering_DX[c] = numbering_dx
                numbering_DY[c] = numbering_dy

            distributions = self._sf_.mesh.elements.distributions

            DisDy = list()
            DisDx = list()

            for core, dis in enumerate(distributions):
                dis_dy = dict()
                dis_dx = dict()
                for c in dis:
                    dis_dy[c] = numbering_DY[c]
                    dis_dx[c] = numbering_DX[c]
                DisDy.append(dis_dy)
                DisDx.append(dis_dx)
        else:
            DisDy = None
            DisDx = None

        DisDx = COMM.scatter(DisDx, root=MASTER_RANK)
        DisDy = COMM.scatter(DisDy, root=MASTER_RANK)

        GV_dict = dict()
        for c in self._sf_.mesh.elements.indices:
            full_vec = np.concatenate([DisDx[c].ravel('F'), DisDy[c].ravel('F')])
            GV_dict[c] = Gathering_Vector(c, full_vec)

        GM = Gathering_Matrix(GV_dict, mesh_type='miUsGrid_TriangularMesh')

        return GM





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
