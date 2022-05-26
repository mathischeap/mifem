# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 2:21 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.base import FrozenOnly

from objects.mpRfT._2d.mesh.rcWds.helpers.vector.rgW import mpRfT2_Mesh_rcWds_Vector_rgW
from objects.mpRfT._2d.mesh.rcWds.helpers.vector.bcW import mpRfT2_Mesh_rcWds_Vector_bcW
from objects.mpRfT._2d.mesh.rcWds.helpers.vector.visualize import mpRfT2_Mesh_rcWds_Vector_Visualize


class mpRfT2_Mesh_rcWds_Vector(FrozenOnly):
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
            The distribution of the data, the classname of the coo_map on which it was built.
        full : bool
            If it covers all root-cells.

        Returns
        -------

        """

        assert isinstance(data, dict), f"data must be in a dict."
        for i in data: assert len(data[i]) == 2

        self._melt_self_()
        self._data_ = data
        self._ndim_ = ndim
        self._distribution_ = distribution
        self._isfull_ = full
        self._freeze_self_()

        return self

    def __getitem__(self, rp):
        return self._data_[rp]

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_Mesh_rcWds_Vector_Visualize(self)
        return self._visualize_

    @property
    def bcW(self):
        return mpRfT2_Mesh_rcWds_Vector_bcW(self)

    @property
    def rgW(self):
        return mpRfT2_Mesh_rcWds_Vector_rgW(self)





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/rcWds/helpers/vector/main.py
    pass
