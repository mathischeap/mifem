# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/24 2:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
import numpy as np

from pyevtk.hl import unstructuredGridToVTK
from root.config.main import RANK, MASTER_RANK, COMM



class miUsTriangular_iS1F_Export_VTK(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.triangle_and_quad(*args, **kwargs)


    def quadratic_triangle(self, filename=None, data_only=False):
        """In each element, we reconstruct with `xi = et = np.array([-1, 0, 1])` and we use the
        values of element boundary to form a `quadratic_triangle (=22)` VTK cell, which means we
        do not make use of the internal reconstruction value (i.e., the value at mapped (0,0).)

        This is an easy and cheap way to visualize the form. However, for higher degree forms, the
        accuracy may be low.

        Parameters
        ----------
        filename
        data_only

        Returns
        -------

        """
        if filename is None:
            filename = self._sf_.name + '_quadratic_triangle'

        connections, offsets, types, global_indices, gathering = \
            self._sf_.mesh.miscellaneous.quadratic_triangle_data

        xi = np.array([-1, 0, 1])
        et = np.array([-0.99, 0, 1])

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
                                      cell_types=types * np.ones_like(offsets),
                                      pointData=pointData
                                      )

    def triangle_and_quad(self, p=None, filename=None, data_only=False):
        """For example, when p = 3, we will have 3*3=9 cells. The three cells close to the singular
        vertex will be triangles(=5), and the rests will be quad(=9).

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
            p = self._sf_.p

        connections, offsets, types, global_indices, gathering = \
            self._sf_.mesh.miscellaneous.triangle_and_quad_data(p)

        xi = np.linspace(-1, 1, p+1)
        et = np.linspace(-0.99, 1, p+1)

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
                return (x, y, z), connections, offsets, types, pointData
            else:
                unstructuredGridToVTK(filename,
                                      x, y, z,
                                      connectivity=connections, offsets=offsets,
                                      cell_types=types,
                                      pointData=pointData
                                      )

if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_1/inner/export/vtk/main.py
    from __init__ import miTri

    fc = miTri.call('st16', 2)

    def p_func(t, x, y): return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) + t
    def q_func(t, x, y): return np.cos(2 * np.pi * x) * np.sin(2 * np.pi * y) + t

    v = fc('vector', [p_func,q_func])

    f1 = fc('1-f-i')
    f1.CF = v
    v.current_time = 0
    f1.discretize()

    f1.export.vtk()

