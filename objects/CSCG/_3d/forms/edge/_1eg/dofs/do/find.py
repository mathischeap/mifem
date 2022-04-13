

from screws.freeze.base import FrozenOnly
from root.config.main import cOmm


class _3dCSCG_E1F_Dofs_DoFind(FrozenOnly):
    """"""
    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._freeze_self_()


    def edge_element_and_local_index_of_dof(self, i):
        """Return the edge_element, and local index (edge-element-wise) of the edge 1-form dof `i`
        in all cores.

        Parameters
        ----------
        i : int
            the dof globally numbered `i` of this 1-edge form.

        Returns
        -------

        """
        EES = self._dofs_._ef_.mesh.edge.elements
        GM = self._dofs_._ef_.numbering.gathering

        mesh_elements_local_numbering = GM.do.find.elements_and_local_indices_of_dof(i)

        if mesh_elements_local_numbering is not None:
            mesh_elements, index = mesh_elements_local_numbering
            mesh_elements = mesh_elements[0]
            index = index[0]

            p = self._dofs_._ef_.space.p

            px, py, pz = p

            if 0 <= index < 4 * px: # N-S edge
                if index < px:
                    edge = "WB"
                    local_index = index
                elif index < 2 * px:
                    edge = "EB"
                    local_index = index - px
                elif index < 3 * px:
                    edge = "WF"
                    local_index = index - 2*px
                else:
                    edge = "EF"
                    local_index = index - 3 * px

            else:
                index -= 4 * px
                if 0 <= index < 4 * py: # W-E edge
                    if index < py:
                        edge = "NB"
                        local_index = index
                    elif index < 2 * py:
                        edge = "SB"
                        local_index = index - py
                    elif index < 3 * py:
                        edge = "NF"
                        local_index = index - 2 * py
                    else:
                        edge = "SF"
                        local_index = index - 3 * py

                else:
                    index -= 4 * py
                    if 0 <= index < 4 * pz: # B-F edge
                        if index < pz:
                            edge = "NW"
                            local_index = index
                        elif index < 2 * pz:
                            edge = "SW"
                            local_index = index - pz
                        elif index < 3 * pz:
                            edge = "NE"
                            local_index = index - 2 * pz
                        else:
                            edge = "SE"
                            local_index = index - 3 * pz
                    else:
                        raise Exception()

            edge_index = ['WB', 'EB', 'WF', 'EF', 'NB', 'SB',
                          'NF', 'SF', 'NW', 'SW', 'NE', 'SE'].index(edge)

            E_MAP = EES.map[mesh_elements]
            edge_element = E_MAP[edge_index]

            edge_element_and_local_index = [edge_element, local_index]
        else:
            edge_element_and_local_index = None

        edge_element_and_local_index = cOmm.allgather(edge_element_and_local_index)

        EE__LI = None
        for ___ in edge_element_and_local_index:
            if ___ is not None:
                if EE__LI is None:
                    EE__LI = ___
                else:
                    assert EE__LI == ___

        return EE__LI