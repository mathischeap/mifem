# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/9/2022 12:31 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from importlib import import_module
from root.config.main import COMM



class miUsGrid_RealTime_TriangularMeshAllocator(FrozenOnly):
    """These meshes are usually structured, and thus we can make them in real time."""

    def __init__(self):
        """"""
        self._freeze_self_()

    @classmethod
    def make_mesh(cls, ID, **kwargs):
        """"""
        path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.' + ID

        function = getattr(import_module(path), ID)

        COMM.barrier()
        boundaries = function(**kwargs)
        source_path = str(__file__).split('allocator')[0] + '__real_time_mesh__.vtu'
        COMM.barrier()

        return source_path, boundaries

    @classmethod
    def check_mesh(cls, ID):
        """Check if the ID refers to a real time triangular mesh."""
        assert ID != '', f"ID can not be empty str."
        if ID in (
            'square', # kwargs: K, d;  a structured mesh in a square [-d/2, d/2] of 2*K**2 elements
        ):
            return True
        else:
            return False






if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/__test__/realtime_meshes/allocator.py
    from __init__ import miTri

    fc = miTri.call('square', 2, K=13, d=3)

    fc.mesh.visualize()