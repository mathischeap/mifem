# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 3:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_MeshIDS_Scalar_Visualize(FrozenOnly):
    """"""

    def __init__(self, scalar):
        """"""
        self._scalar_ = scalar
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    def matplot(self, xy, density=10, title=False):
        """Plot this scalar on the coordinates `xy`.

        Parameters
        ----------
        xy
        density
        title

        Returns
        -------

        """
        assert xy.__class__.__name__ == '_2nCSCG_MRF2_IDS_Vector' and \
               xy.signature == self._scalar_.signature and \
               xy.signature == self._scalar_.mesh.signature, \
            f"coordinates, value and mesh do not match."








if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
