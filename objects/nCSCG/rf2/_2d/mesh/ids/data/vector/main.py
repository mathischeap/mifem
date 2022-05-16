# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 2:21 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.mesh.ids.data.base import _2nCSCG_MRF2_IDS_DataBase
from objects.nCSCG.rf2._2d.mesh.ids.data.vector.visualize import _2nCSCG_MeshIDS_Vector_Visualize
from objects.nCSCG.rf2._2d.mesh.ids.data.vector.BCW import _2nCSCG_MeshIDS_Vector_BCW
from objects.nCSCG.rf2._2d.mesh.ids.data.vector.RGW import _2nCSCG_MeshIDS_Vector_ReGionW


class _2nCSCG_MRF2_IDS_Vector(_2nCSCG_MRF2_IDS_DataBase):
    """"""

    def __init__(self, mesh, data, ndim, distribution, full):
        """"""
        super(_2nCSCG_MRF2_IDS_Vector, self).__init__(mesh, data, ndim, distribution, full)
        for i in self.data: assert len(self.data[i]) == 2
        self._visualize_ = None
        self._freeze_self_()

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2nCSCG_MeshIDS_Vector_Visualize(self)
        return self._visualize_

    @property
    def BCW(self):
        return _2nCSCG_MeshIDS_Vector_BCW(self)

    @property
    def RGW(self):
        return _2nCSCG_MeshIDS_Vector_ReGionW(self)


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
