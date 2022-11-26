# -*- coding: utf-8 -*-
"""
Mainly for visualization reconstruction.

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 4:59 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.mesh.coo_map.helpers.base import mpRfT2_CooMapBase
import numpy as np


class mpRfT2_Mesh_UniformCooMap(mpRfT2_CooMapBase):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._ndim_ = None
        self._1d_nodes = None
        self._2d_nodes = None
        self._coo_ = dict()
        self._freeze_self_()

    @property
    def distribution(self):
        return 'uniform'

    @property
    def ___Pr_rcMC_key___(self):
        return self._mesh_.__Pr_indices_metric_N_key___

    def ___Pr_rcMC_nodes___(self, rp):
        return self[rp]

    @property
    def ___Pr_sgMC_key___(self):
        """A key implying the value for metric involved computing in each root-cell."""
        raise NotImplementedError()

    def ___Pr_sgMC_nodes___(self, rp):
        raise NotImplementedError()

    def __call__(self, factor, ndim=1):
        """

        Parameters
        ----------
        factor : int
            the factor.
        ndim :
            We return `ndim`-dimensional coordinates.

        Returns
        -------

        """
        assert factor > 0 and factor % 1 == 0, f"f={factor} is wrong."
        assert ndim in (1, 2), f"ndim={ndim} is wrong, must be in (1, 2)."
        self._ndim_ = ndim
        self._1d_nodes = np.linspace(-1, 1, factor)
        self._2d_nodes = np.meshgrid(self._1d_nodes, self._1d_nodes, indexing='ij')
        return self

    def __getitem__(self, rp):
        """

        Parameters
        ----------
        rp :
            the __repr__ of a cell.

        Returns
        -------

        """
        #------------ root-cell -----------------------------------------------------------------
        if isinstance(rp, str):

            if '-' not in rp:

                if '' in self._coo_:
                    pass
                else:
                    if self._ndim_ == 1:
                        self._coo_[''] = self._1d_nodes, self._1d_nodes
                    else:
                        self._coo_[''] = self._2d_nodes

                return self._coo_['']

            else:
                ind = rp.split('-')[1]
                if ind in self._coo_:
                    pass
                else:
                    nodes = self._1d_nodes

                    origin, delta = self._mesh_[rp].coordinate_transformation.origin_and_delta

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

        #------------ segment -----------------------------------------------------------------
        elif rp.__class__.__name__ == 'mpRfT2_Segment':
            RP = rp.__repr__()[3]
            if RP == 'c':
                ind = rp.__repr__().split(':')[-1]
            elif RP == 't':
                ind =  rp.__repr__().split('-')[-1]
            else:
                raise Exception()

            if ind in self._coo_:
                pass
            else:
                origin, delta = rp.origin_and_delta

                if isinstance(origin, tuple): # an internal segment
                    direction = ind[0]
                    if direction == 'x':
                        origin = origin[0]
                    elif direction == 'y':
                        origin = origin[1]
                    else:
                        raise Exception()
                else:
                    pass

                end = origin + delta

                nodes = self._1d_nodes
                if origin == -1:
                    _nodes = nodes[nodes <= end]
                else:
                    _nodes = nodes[nodes > origin]
                    _nodes = _nodes[_nodes <= end]

                if len(_nodes) == 0:
                    _nodes = np.array([])
                else:
                    _nodes = 2 * (_nodes - origin) / delta - 1

                self._coo_[ind] = _nodes

            return self._coo_[ind]

        #------------ else -----------------------------------------------------------------
        else:
            raise NotImplementedError(f'Not implemented uniform coo_map for {rp}')
        #===================================================================================



if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/coo_map/helpers/uniform.py
    pass
