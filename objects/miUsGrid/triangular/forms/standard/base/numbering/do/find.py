# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/29/2022 10:42 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from functools import lru_cache


class miUsTriangle_Numbering_DoFind(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def dofs_on_element_edge(self, element, edge):
        """

        Parameters
        ----------
        element : int
            The number of a local element.
        edge : {int, 'D', 'U', 'R'}
            {0, 1, 2}. 0: down edge, 1: upper edge, 2: right edge

        Returns
        -------

        """

        if self._sf_.k == 2:
            raise Exception(f"volume-form does not have dofs on element edge")
        else:
            pass

        local_dofs = self.local_dofs_on_element_edge(edge)
        full_vector = self._sf_.numbering.gathering[element].full_vector

        return full_vector[local_dofs]


    @lru_cache(maxsize=3)
    def local_dofs_on_element_edge(self, edge):
        """

        Parameters
        ----------
        edge : {int, 'D', 'U', 'R'}
            {0, 1, 2}. 0: down edge, 1: upper edge, 2: right edge

        Returns
        -------

        """

        if self._sf_.k == 2:
            raise Exception(f"volume-form does not have dofs on element edge")
        else:
            pass

        if edge == 'D':
            edge = 0
        elif edge == 'U':
            edge = 1
        elif edge == 'R':
            edge = 2
        else:
            assert edge in [0, 1, 2], f"edge={edge} wrong, be in 'DUR' or [0,1,2]."

        if self._sf_.k == 1:
            if self._sf_.__class__.__name__ == 'miUsTriangular_S1F_Outer':
                dy, dx = getattr(self._sf_.space.local_numbering, self._sf_.__class__.__name__)

                if edge == 1: # U edge
                    local_dofs = dy[0, :]
                elif edge == 0: # D edge
                    local_dofs = dy[-1, :]
                elif edge == 2: # edge
                    local_dofs = dx[:, -1]
                else:
                    raise Exception()

                return local_dofs

            elif self._sf_.__class__.__name__ == 'miUsTriangular_S1F_Inner':
                dx, dy = getattr(self._sf_.space.local_numbering, self._sf_.__class__.__name__)

                if edge == 1: # U edge
                    local_dofs = dy[0, :]
                elif edge == 0: # D edge
                    local_dofs = dy[-1, :]
                elif edge == 2: # edge
                    local_dofs = dx[:, -1]
                else:
                    raise Exception()

                return local_dofs

            else:
                raise NotImplementedError()

        elif self._sf_.k == 0:
            Local_numbering = getattr(self._sf_.space.local_numbering, self._sf_.__class__.__name__)[0]

            if edge == 1:
                return Local_numbering[0, :]
            elif edge == 0:
                return Local_numbering[-1, :]
            elif edge == 2:
                return Local_numbering[:, -1]
            else:
                raise Exception()

        else:
            raise Exception()





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
