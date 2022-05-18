# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 7:11 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np

from objects.nCSCG.rf2._2d.mesh.coordinates.distributions.base import _2nCSCG_MRF2_CooDistributionBase


class Homogeneous(_2nCSCG_MRF2_CooDistributionBase):
    """"""

    def __init__(self, mesh):
        """"""
        super(Homogeneous, self).__init__(mesh)
        self._distribution_ = 'homogeneous'
        self._1d_nodes = None
        self._2d_nodes = None
        self._coo_ = dict()
        self._freeze_self_()

    def __call__(self, f, ndim=1):
        """

        Parameters
        ----------
        f : int
            the factor.
        ndim :
            We return `ndim`-dimensional coordinates.

        Returns
        -------

        """
        assert self.signature == self._mesh_.signature, f"signature dis-match"
        assert f > 0 and f % 1 == 0, f"f={f} is wrong."
        assert ndim in (1, 2), f"ndim={ndim} is wrong, must be in (1, 2)."
        self._ndim_ = ndim

        self._1d_nodes = np.linspace(-1, 1, f)
        self._2d_nodes = np.meshgrid(self._1d_nodes, self._1d_nodes, indexing='ij')

        return self

    def __getitem__(self, indices):
        """

        Parameters
        ----------
        indices :
            the indices of a cell.

        Returns
        -------

        """

        level = len(indices)
        if level == 1:

            if 0 in self._coo_:
                pass
            else:
                if self._ndim_ == 1:
                    self._coo_[0] = self._1d_nodes, self._1d_nodes
                else:
                    self._coo_[0] = self._2d_nodes

            return self._coo_[0]

        else:
            ind = indices[1:]
            if ind in self._coo_:
                pass
            else:
                nodes = self._1d_nodes

                origin, delta = self._mesh_.do.find.origin_and_delta(indices)

                ox, oy = origin
                ex, ey = ox + delta, oy + delta

                if ox == -1:
                    x_nodes = nodes[nodes <= ex]
                else:
                    x_nodes = nodes[nodes > ox]
                    x_nodes = x_nodes[x_nodes <= ex]

                if oy == -1:
                    y_nodes = nodes[nodes <= ey]
                else:
                    y_nodes = nodes[nodes > oy]
                    y_nodes = y_nodes[y_nodes <= ey]

                if len(x_nodes) == 0:
                    x_nodes = np.array([])
                else:
                    x_nodes = 2 * (x_nodes - ox) / delta - 1

                if len(y_nodes) == 0:
                    y_nodes = np.array([])
                else:
                    y_nodes = 2 * (y_nodes - oy) / delta - 1

                if self._ndim_ == 1:
                    self._coo_[ind] = x_nodes, y_nodes
                else:
                    self._coo_[ind] = np.meshgrid(x_nodes, y_nodes, indexing='ij')

            return self._coo_[ind]










if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass