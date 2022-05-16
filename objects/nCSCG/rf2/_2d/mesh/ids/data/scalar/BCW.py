# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 4:23 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np

class _2nCSCG_MeshIDS_Scalar_BCW(FrozenOnly):
    """Data are arranged Base-Cell-Wise"""

    def __init__(self, scalar):
        """"""
        # below is an over-powerful constraint, we only need base-cell-wise homogeneous, but below is globally homogeneous.
        assert scalar.distribution == 'homogeneous', \
            f"distribution={scalar.distribution} wrong, only accept homogeneous distribution."
        assert scalar.ndim == 2, f"Only 2-dim data can be BCW arranged."
        self._scalar_ = scalar
        assert scalar.signature == scalar.mesh.signature, f"mesh has been changed."
        self.___Pr_arrange_data___()
        self._freeze_self_()

    def ___Pr_arrange_data___(self):
        """"""
        mesh = self._scalar_.mesh
        self._BCW_ = dict()
        for i in mesh.base_cells:
            self._BCW_[i] = self.___Pr_find_data_for_cell___(i)

    def ___Pr_find_data_for_cell___(self, i):
        """"""
        cell = self._scalar_.mesh(i)
        data = self._scalar_.data
        if cell.IS.root:
            try:
                return data[cell.__repr__()][0]
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
                return np.bmat([(d0, d2),
                                (d1, d3)])

    def __getitem__(self, i):
        """

        Parameters
        ----------
        i : int
            Local base-cell number.

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
