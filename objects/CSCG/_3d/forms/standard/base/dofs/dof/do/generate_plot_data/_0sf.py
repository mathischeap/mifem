# -*- coding: utf-8 -*-
import numpy as np

from screws.freeze.base import FrozenOnly
from root.config.main import COMM


class GPD_0SF(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._freeze_self_()

    def __call__(self, mesh_element, zoom=1):
        """"""
        sf = self._dof_._dofs_._sf_
        mesh = sf.mesh

        #------- parse the mesh elements to be checked -----------------------------------


        if mesh_element is None: # we will return the data from all possible local mesh elements.
            MESH_ELEMENTS = list()

            positions = self._dof_.positions

            for pos in positions:
                mesh_element, _ = pos
                assert mesh_element in mesh.elements

                MESH_ELEMENTS.append(mesh_element)

            MESH_ELEMENTS = COMM.allgather(MESH_ELEMENTS)
            ___ = set()
            for _ in MESH_ELEMENTS:
                ___.update(_)
            MESH_ELEMENTS = list(___)

        elif isinstance(mesh_element, int):
            MESH_ELEMENTS = [mesh_element,]

        else:
            raise Exception()

        #---------------------------------------------------------------------------------------
        LOCAL_POSITIONS = list()
        positions = self._dof_.positions
        for pos in positions:
            LOCAL_POSITIONS.append(pos[0])
        nx, ny, nz = sf.space.nodes

        local_numbering = sf.numbering.local[0]
        COO = dict()
        for me in MESH_ELEMENTS:

            if me in LOCAL_POSITIONS:

                COO[me] = list()

                for pos in positions:
                    if me == pos[0]:

                        local_index = pos[1]

                        i, j, k = np.argwhere(local_numbering == local_index)[0]

                        x = nx[i] * zoom
                        y = ny[j] * zoom
                        z = nz[k] * zoom

                        mesh_element = mesh.elements[me]

                        x, y, z = mesh_element.coordinate_transformation.mapping(x,y,z)

                        COO[me].extend([x[0], y[0], z[0]])

                        break


        RETURN = dict()
        RETURN['DOF_COO'] = COO # coordinates

        return RETURN