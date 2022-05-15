# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 2:21 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.mesh.ids.data.base import _2nCSCG_MRF2_IDS_DataBase


class _2nCSCG_MRF2_IDS_Vector(_2nCSCG_MRF2_IDS_DataBase):
    """"""

    def __init__(self, mesh, data, distribution, full):
        """"""
        super(_2nCSCG_MRF2_IDS_Vector, self).__init__(mesh, data, distribution, full)
        self._freeze_self_()





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
