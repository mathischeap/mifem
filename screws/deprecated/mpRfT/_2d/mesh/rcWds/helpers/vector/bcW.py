# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 6:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np

class mpRfT2_Mesh_rcWds_Vector_bcW(FrozenOnly):
    """Data are arranged Base-Cell-Wise"""

    def __init__(self, vector):
        """"""
        # below is an over-powerful constraint, we only need base-cell-wise uniform, but below is globally uniform.
        assert vector._distribution_ == 'uniform', \
            f"distribution={vector._distribution_} wrong, only accept uniform distribution."
        assert vector._ndim_ == 2, f"Only 2-dim data can be basic-cell-wise arranged."
        assert vector._isfull_, f"only full rcWds vector can be reshaped into bcW!"
        self._vector_ = vector
        self.___Pr_arrange_data___()
        self._freeze_self_()

    def ___Pr_arrange_data___(self):
        """"""
        mesh = self._vector_._mesh_
        self._BCW_ = dict()
        for i in mesh.basic_cells:
            self._BCW_[i] = self.___Pr_find_data_for_cell___(i)

    def ___Pr_find_data_for_cell___(self, i):
        """"""
        cell = self._vector_._mesh_[i]
        data = self._vector_._data_
        if cell.IS.root:
            try:
                return data[cell.__repr__()]
            except KeyError:
                return None
        else:
            sc0, sc1, sc2, sc3 = [cell.indices + (_,) for _ in range(4)]
            d0 = self.___Pr_find_data_for_cell___(sc0)
            d1 = self.___Pr_find_data_for_cell___(sc1)
            d2 = self.___Pr_find_data_for_cell___(sc2)
            d3 = self.___Pr_find_data_for_cell___(sc3)
            if any([d is None for d in [d0, d1, d2, d3]]):
                return None
            else:
                d00, d01 = d0
                d10, d11 = d1
                d20, d21 = d2
                d30, d31 = d3
                return np.bmat([(d00, d20), (d10, d30)]), np.bmat([(d01, d21), (d11, d31)])

    def __getitem__(self, i):
        """

        Parameters
        ----------
        i : int
            Local basic-cell number.

        Returns
        -------

        """
        return self._BCW_[i]

    def __iter__(self):
        for i in self._BCW_:
            yield i


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
