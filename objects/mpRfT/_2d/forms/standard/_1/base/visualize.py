# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/24/2022 8:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_S1F_Visualize(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.matplot(*args, **kwargs)

    def matplot(self, density=20, plot_type='contourf', **kwargs):
        """"""
        mesh = self._f_.mesh
        coo_map = mesh.coo_map.uniform(density, ndim=1)
        xy, v = self._f_.reconstruct(coo_map, ravel=False)
        v.visualize(xy, plot_type=plot_type, **kwargs)



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
