# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/13/2022 10:44 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
import vtk
# noinspection PyUnresolvedReferences
from vtk.util.numpy_support import vtk_to_numpy
from root.config.main import cOmm, mAster_rank, rAnk, sIze, np
from screws.functions._2d_space.angle import angle



class miUsGrid_TriangularMesh_Construct(FrozenOnly):
    """
    A triangle cell has three vertexes, there are number 0, 1, 2, 0->1->2->0 must form a right hand
    rule.

    The three edges of a cell are: 0->1, 0->2, 1->2. This is very important.

    """

    def __init__(self, source):
        """

        Parameters
        ----------
        source
        """

        if isinstance(source, str) and source[-4:] == '.vtu': # unstructured VTK file.
            self._cells_and_points_ = self.from_vtu_file(source)
        else:
            raise NotImplementedError(f"")

        self._freeze_self_()

    def __call__(self, mesh):
        """We analyze the cells and points data."""
        if rAnk == mAster_rank:
            cells, points = self._cells_and_points_
            NUM_element = cells.shape[0]

            #--- we go through the cells to check if we need to switch vertex 1 and vertex 2.------1
            for i, cell in enumerate(cells):
                c0, c1, c2 = points[cell]

                angle0 = angle(c0, c1)
                angle1 = angle(c0, c2)

                if angle0 == 0:
                    if 0 < angle1 < np.pi:
                        pass
                    elif angle1 > np.pi:# circle 0->1->2->0 is left-hand-ruled.
                        # switch vertex 1 and vertex 2
                        cells[i][1:] = cells[i][1:][::-1]

                    else:
                        raise Exception(f"cell[{i}] is illegal.")

                elif angle1 == 0:
                    if angle0 > np.pi:
                        pass
                    elif 0 < angle0 < np.pi:# circle 0->1->2->0 is left-hand-ruled.
                        # switch vertex 1 and vertex 2
                        cells[i][1:] = cells[i][1:][::-1]
                    else:
                        raise Exception(f"cell[{i}] is illegal.")

                else:
                    ad = angle0 - angle1
                    if ad == 0:
                        raise Exception(f"cell #[{i}] is not a triangle.")
                    elif ad > 0 : # circle 0->1->2->0 is left-hand-ruled.
                        # switch vertex 1 and vertex 2
                        cells[i][1:] = cells[i][1:][::-1]

                    else:
                        pass

            #--------- now make the element map in the master core --------------------------------1
            ELEMENT_MAP = dict()
            for i in range(NUM_element):
                ELEMENT_MAP[i] = ['', '', '']

            for i, cell in enumerate(cells):
                v0, v1, v2 = cell

                edges = [(v0, v1),
                         (v0, v2),
                         (v1, v2)]

                for j in range(i+1, NUM_element):
                    target_vertexes = cells[j]
                    num_shared_vertexes = 0
                    for v in cell:
                        if v in target_vertexes:
                            num_shared_vertexes += 1

                    if num_shared_vertexes == 2:
                        tv0, tv1, tv2 = target_vertexes
                        Tes = {0:(tv0, tv1),
                                1:(tv0, tv2),
                                2:(tv1, tv2),
                                -1:(tv1, tv0),
                                -2:(tv2, tv0),
                                -3:(tv2, tv1)}

                        for m, E in enumerate(edges):
                            for n in Tes:
                                Te = Tes[n]
                                if E == Te:
                                    break
                            # noinspection PyUnboundLocalVariable
                            if E == Te:
                                break

                        # noinspection PyUnboundLocalVariable
                        if n == -1:
                            n = 0
                            sign = '-'
                        elif n == -2:
                            n = 1
                            sign = '-'
                        elif n == -3:
                            n = 2
                            sign = '-'
                        else:
                            sign = '+'

                        mark_for_i = str(j) + sign + str(n)
                        # noinspection PyUnboundLocalVariable
                        mark_for_j = str(i) + sign + str(m)

                        ELEMENT_MAP[i][m] = mark_for_i
                        ELEMENT_MAP[j][n] = mark_for_j

                    elif num_shared_vertexes <= 1:
                        pass

                    elif '' not in ELEMENT_MAP[i]:
                        break
                    else:
                        raise Exception()

            #------------------- distribute the cells to cores ------------------------------------1
            distribution = [NUM_element // sIze + (1 if x < NUM_element % sIze else 0)
                            for x in range(sIze)][::-1] # `i`th core will handle the amount of distribution[i] elements.
            element_maps_DIS = [dict()  for _ in range(sIze)]
            LOCAL_ELEMENT_RANGES = list()
            CELLS = [dict()  for _ in range(sIze)]
            POINTS = [dict()  for _ in range(sIze)]

            num_points = len(points)

            current_num = 0
            for i, num in enumerate(distribution):
                RANGE = range(current_num, current_num + num)
                for c in RANGE:
                    element_maps_DIS[i][c] = ELEMENT_MAP[c]

                    CELLS[i][c] = list(cells[c])
                    for cell in cells[c]:
                        POINTS[i][cell] = list(points[cell])

                current_num += num
                LOCAL_ELEMENT_RANGES.append(RANGE)

        else:
            element_maps_DIS = None
            LOCAL_ELEMENT_RANGES = None
            CELLS = None
            POINTS = None
            num_points = 0

        CELLS = cOmm.scatter(CELLS, root=mAster_rank)
        POINTS = cOmm.scatter(POINTS, root=mAster_rank)
        element_maps = cOmm.scatter(element_maps_DIS, root=mAster_rank)
        distribution = cOmm.bcast(LOCAL_ELEMENT_RANGES, root=mAster_rank)
        LOCAL_ELEMENT_RANGES = distribution[rAnk]

        num_points = cOmm.bcast(num_points, root=mAster_rank)

        mesh._elements_._cells_ = CELLS
        mesh._elements_._points_ = POINTS
        mesh._elements_._map_ = element_maps
        mesh._elements_._distributions_ = distribution
        mesh._elements_._range_ = LOCAL_ELEMENT_RANGES
        mesh._elements_.__num_GLOBAL_points__ = num_points



    @staticmethod
    def from_vtu_file(source_file_name):
        """"""
        if rAnk == mAster_rank:
            # noinspection PyUnresolvedReferences
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(source_file_name)
            reader.Update()
            data = reader.GetOutput()
            points = vtk_to_numpy(data.GetPoints().GetData())[:,:2]
            cells = vtk_to_numpy(data.GetCells().GetData())
            number_of_cells  = data.GetNumberOfCells()
            cell_types = cells[::4]
            assert np.all(cell_types == 3), f"Not a triangular mesh!"

            cells = cells.reshape((number_of_cells, 4), order='C')[:, 1:]
            return cells, points
        else:
            return None



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/construct/main.py
    pass
