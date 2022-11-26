# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 7:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_ScalarVisualize(FrozenOnly):
    """"""

    def __init__(self, cf):
        """"""
        self._cf_ = cf
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.matplot(*args, **kwargs)


    def matplot(self, *args, **kwargs):
        """"""
        ftype = self._cf_.ftype

        if ftype == 'standard':
            return self.___Pr_matplot_standard___(*args, **kwargs)
        else:
            raise NotImplementedError(f"{ftype}.")


    def ___Pr_matplot_standard___(self, density=20, plot_type='contourf', **kwargs):
        """"""
        mesh = self._cf_.mesh
        coo = mesh.coo_map.uniform(density, ndim=2)
        xy, v = self._cf_.reconstruction(coo, ravel=False)
        v.visualization(xy, plot_type=plot_type, **kwargs)







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/cf/scalar/visualize.py
    pass