

from screws.freeze.base import FrozenOnly
import numpy as np
from root.config.main import cOmm

class _3dCSCG_1EF_Special(FrozenOnly):
    """"""
    def __init__(self, ef):
        self._ef_ = ef
        self._freeze_self_()


    def generate_plot_data_for_dof(self, i, density=5, zoom=1):
        """Return the plot data in all cores.

        Parameters
        ----------
        i
        density
        zoom : float
            It will only apply to the direction of the edge element (where the edge dof is)
            (one of N-S, W-E, or B-F).

        Returns
        -------

        """
        GM = self._ef_.numbering.gathering

        elements_indices = GM.do.find.elements_and_local_indices_of_dof(i)

        if elements_indices is not None: # find the data in all cores that involve this edge dof.
            elements, indices = elements_indices

            element = elements[0]
            index = indices[0]

            ME = self._ef_.mesh.elements[element]

            p = self._ef_.space.p
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
                        raise Exception

            P1 =  np.ones(density)
            M1 = - np.ones(density)

            nodes = self._ef_.space.nodes

            if edge == "WB":
                node = nodes[0]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = N, M1, M1
            elif edge == "EB":
                node = nodes[0]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = N, P1, M1
            elif edge == "WF":
                node = nodes[0]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = N, M1, P1
            elif edge == "EF":
                node = nodes[0]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = N, P1, P1


            elif edge == 'NB':
                node = nodes[1]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = M1, N, M1
            elif edge == 'SB':
                node = nodes[1]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = P1, N, M1
            elif edge == 'NF':
                node = nodes[1]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = M1, N, P1
            elif edge == 'SF':
                node = nodes[1]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = P1, N, P1

            elif edge == 'NW':
                node = nodes[2]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = M1, M1, N
            elif edge == 'SW':
                node = nodes[2]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = P1, M1, N
            elif edge == 'NE':
                node = nodes[2]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = M1, P1, N
            elif edge == 'SE':
                node = nodes[2]
                N = np.linspace(node[local_index], node[local_index+1], density) * zoom
                x, y, z = P1, P1, N

            else:
                raise Exception()

            xyz = ME.coordinate_transformation.mapping(x, y, z)

        else:
            xyz = None

        XYZ = cOmm.allgather(xyz)
        xyz = list()
        for ___ in XYZ:
            if ___ is not None:
                xyz.append(___)

        return xyz