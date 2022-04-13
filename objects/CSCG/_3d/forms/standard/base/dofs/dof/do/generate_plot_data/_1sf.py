

from screws.freeze.base import FrozenOnly
import numpy as np
from root.config.main import cOmm, MPI


class GPD_1SF(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._freeze_self_()

    def __call__(self, mesh_element, density=5, zoom=1):
        """We generate the plotting data considering the dof is in the
        mesh element #`'`mesh_element`.

        Parameters
        ----------
        mesh_element
        density

        Returns
        -------
        RETURN : dict

            'EFD' : Element Frame Data
        """
        FIND = 0
        positions = self._dof_.positions
        for pos in positions:
            _ = pos[0]
            if _ == mesh_element:
                FIND += 1

        FIND = cOmm.allreduce(FIND, op=MPI.SUM)
        assert FIND == 1, f"I (dof # {self._dof_.i}) not in the mesh element #{mesh_element}"

        if mesh_element in self._dof_._sf_.mesh.elements:
            EFD = self._dof_._sf_.___PRIVATE_element_grid_data_generator_1___(
                mesh_element,
                density=2*density,
                zoom=zoom)

            for pos in positions:
                _ = pos[0]
                if _ == mesh_element:
                    index = pos[1]
                    I, j, k = np.where(self._dof_._sf_.numbering.local[0] == index)
                    if len(I) == 0:
                        I, j, k = np.where(self._dof_._sf_.numbering.local[1] == index)
                        if len(I) == 0:
                            I, j, k = np.where(self._dof_._sf_.numbering.local[2] == index)
                            assert len(I) != 0, f"Must have found the index."
                            IND = 2
                        else:
                            IND = 1
                    else:
                        IND = 0
                    I, j, k = I[0], j[0], k[0]
                    nodes = self._dof_._sf_.space.nodes
                    if IND == 0:
                        x, y, z = np.linspace(nodes[0][I], nodes[0][I + 1], density), \
                                  nodes[1][j] * np.ones(density), \
                                  nodes[2][k] * np.ones(density)
                    elif IND == 1:
                        x, y, z = nodes[0][I] * np.ones(density), \
                                  np.linspace(nodes[1][j], nodes[1][j + 1], density), \
                                  nodes[2][k] * np.ones(density)
                    elif IND == 2:
                        x, y, z = nodes[0][I] * np.ones(density), \
                                  nodes[1][j] * np.ones(density), \
                                  np.linspace(nodes[2][k], nodes[2][k + 1], density)
                    else:
                        raise Exception()
                    x *= zoom
                    y *= zoom
                    z *= zoom
                    xyz = self._dof_._sf_._mesh_.elements[mesh_element].\
                              coordinate_transformation.mapping(x, y, z)
                else:
                    pass
        else:
            EFD = None
            xyz = None

        RETURN = dict()
        RETURN['EFD'] = EFD
        RETURN['DOF_COO'] = xyz

        return RETURN
