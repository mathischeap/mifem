# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/14 6:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')
from components.freeze.base import FrozenOnly

from objects.mpRfT._2d.mesh.rcWds.helpers.scalar.rgW import mpRfT2_Mesh_rcWds_Scalar_rgW
from objects.mpRfT._2d.mesh.rcWds.helpers.scalar.bcW import mpRfT2_Mesh_rcWds_Scalar_bcW
from objects.mpRfT._2d.mesh.rcWds.helpers.scalar.visualize import mpRfT2_Mesh_rcWds_Scalar_Visualize

import numpy as np


class mpRfT2_Mesh_rcWds_Scalar(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._visualize_ = None
        self._freeze_self_()

    def __call__(self, data, ndim, distribution, full):
        """

        Parameters
        ----------
        data
        ndim : the dimension of the root-cell-wise mesh.
        distribution :
            The distribution of the data, the classname of the coo_map on which the data were built.
        full : bool
            If it covers all root-cells.

        Returns
        -------

        """
        assert isinstance(data, dict), f"data must be in a dict."
        for i in data: assert len(data[i]) == 1

        self._melt_self_()
        self._data_ = data
        self._ndim_ = ndim
        self._distribution_ = distribution
        self._isfull_ = full
        self._freeze_self_()

        return self

    def __getitem__(self, rp):
        return self._data_[rp]

    def __iter__(self):
        for rp in self._data_:
            yield rp

    @property
    def visualization(self):
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_Mesh_rcWds_Scalar_Visualize(self)
        return self._visualize_

    @property
    def bcW(self):
        return mpRfT2_Mesh_rcWds_Scalar_bcW(self)

    @property
    def rgW(self):
        return mpRfT2_Mesh_rcWds_Scalar_rgW(self)


    def __add__(self, other):
        """"""
        assert other._mesh_ is self._mesh_
        assert other.__class__.__name__ == 'mpRfT2_Mesh_rcWds_Scalar'
        assert other._distribution_ == self._distribution_
        assert other._ndim_ == self._ndim_
        assert other._isfull_ == self._isfull_
        data = dict()
        for i in self._data_:
            data[i] = [self._data_[i][0] + other._data_[i][0], ]
        return mpRfT2_Mesh_rcWds_Scalar(self._mesh_)(data, self._ndim_, self._distribution_, self._isfull_)

    def __sub__(self, other):
        """"""
        assert other._mesh_ is self._mesh_
        assert other.__class__.__name__ == 'mpRfT2_Mesh_rcWds_Scalar'
        assert other._distribution_ == self._distribution_
        assert other._ndim_ == self._ndim_
        assert other._isfull_ == self._isfull_
        data = dict()
        for i in self._data_:
            data[i] = [self._data_[i][0] - other._data_[i][0], ]
        return mpRfT2_Mesh_rcWds_Scalar(self._mesh_)(data, self._ndim_, self._distribution_, self._isfull_)

    def __abs__(self):
        """"""
        data = dict()
        for i in self._data_:
            data[i] = [np.abs(self._data_[i][0]), ]
        return mpRfT2_Mesh_rcWds_Scalar(self._mesh_)(data, self._ndim_, self._distribution_, self._isfull_)


    def pow(self, n):
        data = dict()
        for i in self._data_:
            data[i] = [self._data_[i][0]**n, ]
        return mpRfT2_Mesh_rcWds_Scalar(self._mesh_)(data, self._ndim_, self._distribution_, self._isfull_)

    @property
    def magnitude(self):
        data = dict()
        for i in self._data_:
            di = np.abs(self._data_[i][0])
            di[ di < 1e-10] = 1e-10
            data[i] = [np.log10(di), ]
        return mpRfT2_Mesh_rcWds_Scalar(self._mesh_)(data, self._ndim_, self._distribution_, self._isfull_)



if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/its/data_tree.py
    pass
