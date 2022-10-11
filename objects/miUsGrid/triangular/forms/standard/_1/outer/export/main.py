# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 4:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
import numpy as np

from pyevtk.hl import unstructuredGridToVTK
from root.config.main import rAnk, mAster_rank, cOmm


class miUsTriangular_oS1F_Export(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.vtk_quadratic_triangle(*args, **kwargs)

    def vtk_quadratic_triangle(self, filename=None, data_only=False):
        """"""
        if filename is None:
            filename = self._sf_.name

        connections, offsets, types, global_indices, gathering = \
            self._sf_.mesh.miscellaneous.quadratic_triangle_data

        xi = np.array([-1, 0 ,1])
        et = np.array([-0.99, 0, 1])

        xy, Rc = self._sf_.reconstruct(xi, et, ravel=True)

        xy = cOmm.gather(xy, root=mAster_rank)
        Rc = cOmm.gather(Rc, root=mAster_rank)

        if rAnk == mAster_rank:
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
            vx = np.zeros(num_data+1)
            vy = np.zeros(num_data+1)

            for e in xy:
                _x, _y = xy[e]
                x[gathering[e,:]] = _x[global_indices]
                y[gathering[e,:]] = _y[global_indices]
                _rc = Rc[e]
                vx[gathering[e,:]] = _rc[0][global_indices]
                vy[gathering[e,:]] = _rc[1][global_indices]

            pointData = {
                f'{self._sf_.name}': (vx, vy, z)
            }

            if data_only:
                return (x, y, z), connections, pointData
            else:
                unstructuredGridToVTK(filename,
                                      x, y, z,
                                      connectivity=connections, offsets=offsets,
                                      cell_types=22 * np.ones_like(offsets),
                                      pointData=pointData
                                      )



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_1/outer/export/main.py
    from __init__ import miTri

    fc = miTri.form('st16', 2)

    def p_func(t, x, y): return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) + t
    def q_func(t, x, y): return np.cos(2 * np.pi * x) * np.sin(2 * np.pi * y) + t

    v = fc('vector', [p_func,q_func])

    f1 = fc('1-f-o')
    f1.CF = v
    v.current_time = 0
    f1.discretize()

    f1.export()


