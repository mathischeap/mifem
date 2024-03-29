# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/13/2022 8:53 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
import numpy as np

from pyevtk.hl import unstructuredGridToVTK
from root.config.main import RANK, MASTER_RANK, COMM


class miUsGrid_Triangular_Scalar_Export_VTK(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        return self.triangle_and_quad(*args, **kwargs)

    def triangle_and_quad(self, p=None, filename=None, data_only=False):
        """

        Parameters
        ----------
        p
        filename
        data_only

        Returns
        -------

        """
        if filename is None:
            filename = self._sf_.name + '_triangle_and_quad'
        if p is None:
            p = 3

        connections, offsets, types, global_indices, gathering = \
            self._sf_.mesh.miscellaneous.triangle_and_quad_data(p)

        xi = et = np.linspace(-1, 1, p+1)

        xi, et = np.meshgrid(xi, et, indexing='ij')

        xy, Rc = self._sf_.reconstruct(xi, et, ravel=True)

        xy = COMM.gather(xy, root=MASTER_RANK)
        Rc = COMM.gather(Rc, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            _ = dict()
            for _xy in xy:
                _.update(_xy)
            xy = _
            _ = dict()
            for _Rc in Rc:
                _.update(_Rc)
            Rc = _

            num_data = np.max(connections)
            x = np.zeros(num_data+1)
            y = np.zeros(num_data+1)
            z = np.zeros(num_data+1)
            v = np.zeros(num_data+1)

            for e in xy:
                _x, _y = xy[e]
                x[gathering[e,:]] = _x[global_indices]
                y[gathering[e,:]] = _y[global_indices]
                _rc = Rc[e][0]
                v[gathering[e,:]] = _rc[global_indices]

            pointData = {
                f'{self._sf_.name}': v
            }
            if data_only:
                return (x, y, z), connections, offsets, types, pointData
            else:
                unstructuredGridToVTK(filename,
                                      x, y, z,
                                      connectivity=connections, offsets=offsets,
                                      cell_types=types,
                                      pointData=pointData
                                      )


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
