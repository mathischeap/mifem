# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/14 5:49 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.ids.data.scalar.main import _2nCSCG_MRF2_IDS_Scalar
from objects.nCSCG.rf2._2d.mesh.ids.data.vector.main import _2nCSCG_MRF2_IDS_Vector





class _2nCSCG_MeshRF2_IndicesDataStorage(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, data_type, data, ndim, distribution, full):
        """

        Parameters
        ----------
        data_type
        data
        ndim : int
            The dimensions of the data
        distribution : str
        full : bool
            Whether the data cover all current local (sub-)cells.

        Returns
        -------

        """
        if data_type == 'scalar':
            return _2nCSCG_MRF2_IDS_Scalar(self._mesh_, data, ndim, distribution, full)
        elif data_type == 'vector':
            return _2nCSCG_MRF2_IDS_Vector(self._mesh_, data, ndim, distribution, full)
        else:
            raise NotImplementedError(f"not implemented for data type= {data_type}.")





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
