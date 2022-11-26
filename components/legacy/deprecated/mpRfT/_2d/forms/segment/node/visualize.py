# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/01 11:21 AM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_NSgF_Visualize(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    def matplot(self, dp=3, **kwargs):
        """"""
        mesh = self._t_.mesh
        coo = mesh.coo_map.Lobatto(dp)
        xy, v = self._t_.reconstruction(coo)
        v.visualization(xy, **kwargs)








if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/node/visualize.py
    pass
