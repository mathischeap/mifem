# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/14 6:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.mesh.ids.data.base import _2nCSCG_MRF2_IDS_DataBase
from objects.nCSCG.rf2._2d.mesh.ids.data.scalar.visualize import _2nCSCG_MeshIDS_Scalar_Visualize
from objects.nCSCG.rf2._2d.mesh.ids.data.scalar.do import _2nCSCG_MeshIDS_Scalar_Do
from objects.nCSCG.rf2._2d.mesh.ids.data.scalar.BCW import _2nCSCG_MeshIDS_Scalar_BCW
from objects.nCSCG.rf2._2d.mesh.ids.data.scalar.RGW import _2nCSCG_MeshIDS_Scalar_ReGionW




class _2nCSCG_MRF2_IDS_Scalar(_2nCSCG_MRF2_IDS_DataBase):
    """"""

    def __init__(self, mesh, data, ndim, distribution, full):
        """"""
        super(_2nCSCG_MRF2_IDS_Scalar, self).__init__(mesh, data, ndim, distribution, full)
        for i in self.data: assert len(self.data[i]) == 1
        self._visualize_ = None
        self._do_ = None
        self._freeze_self_()

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2nCSCG_MeshIDS_Scalar_Visualize(self)
        return self._visualize_

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _2nCSCG_MeshIDS_Scalar_Do(self)
        return self._do_

    @property
    def BCW(self):
        return _2nCSCG_MeshIDS_Scalar_BCW(self)

    @property
    def RGW(self):
        return _2nCSCG_MeshIDS_Scalar_ReGionW(self)




if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/its/data_tree.py
    pass
