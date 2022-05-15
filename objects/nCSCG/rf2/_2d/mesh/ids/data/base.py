# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 2:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_MRF2_IDS_DataBase(FrozenOnly):
    """"""

    def __init__(self, mesh, data, distribution, full):
        """

        Parameters
        ----------
        mesh
        data
        distribution
        full : bool
            If the data cover all current local (sub-)cells.
        """
        self._mesh_ = mesh
        self._signature_ = mesh._signature_
        self._distribution_ = distribution
        # this signature implies the mesh on which this tree is built on.
        self._data_ = data
        self._full_ = full

    @property
    def mesh(self):
        return self._mesh_

    @property
    def data(self):
        return self._data_

    @property
    def distribution(self):
        return self._distribution_

    @property
    def signature(self):
        return self._signature_

    @property
    def _isfull_(self):
        return self._full_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
