# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/15 8:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_S1F_Numbering_DO_FIND(FrozenOnly):
    """"""

    def __init__(self, numbering):
        """"""
        self._numbering_ = numbering
        self._f_ = numbering._f_
        self._edge_names_ = ['U', 'D', 'L', 'R']
        self._cache0_ = dict()
        self._freeze_self_()


    def local_numbering_of_dofs_on(self, rc_rp, edge_name):
        """

        Parameters
        ----------
        rc_rp: str
            The repr of a root-cell
        edge_name: str
            {'U', 'D', 'L', 'R'}.

        Returns
        -------

        """
        assert edge_name in self._edge_names_, f"edge_name={edge_name} is wrong."
        N = self._f_.N[rc_rp]

        key = str(N)+edge_name

        if key in self._cache0_:
            pass

        else:
            local_numbering = self._numbering_.local[rc_rp]

            if self._f_.orientation == 'outer':
                if edge_name == 'U':
                    dofs = local_numbering[0][0, :]
                elif edge_name == 'D':
                    dofs = local_numbering[0][-1, :]
                elif edge_name == 'L':
                    dofs = local_numbering[1][:, 0]
                elif edge_name == 'R':
                    dofs = local_numbering[1][:, -1]
                else:
                    raise Exception()

            elif self._f_.orientation == 'inner':
                if edge_name == 'U':
                    dofs = local_numbering[1][0, :]
                elif edge_name == 'D':
                    dofs = local_numbering[1][-1, :]
                elif edge_name == 'L':
                    dofs = local_numbering[0][:, 0]
                elif edge_name == 'R':
                    dofs = local_numbering[0][:, -1]
                else:
                    raise Exception()

            else:
                raise Exception()

            self._cache0_[key]  = dofs

        return self._cache0_[key]







if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
