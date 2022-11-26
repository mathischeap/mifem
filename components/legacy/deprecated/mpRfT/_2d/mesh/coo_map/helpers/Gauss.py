# -*- coding: utf-8 -*-
"""
Mainly for integral reconstruction.

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 5:05 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.mesh.coo_map.helpers.base import mpRfT2_CooMapBase
import numpy as np
from components.quadrature import Quadrature


class mpRfT2_Mesh_GaussCooMap(mpRfT2_CooMapBase):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._dp_ = None
        self._coo_ = dict()
        self._freeze_self_()

    @property
    def distribution(self):
        return 'Gauss'

    @property
    def ___Pr_rcMC_key___(self):
        return self._mesh_.___Pr_metric_N_key___

    def ___Pr_rcMC_nodes___(self, rp):
        return self[rp][0][1]

    @property
    def ___Pr_sgMC_key___(self):
        """A key implying the value for metric involved computing in each root-cell."""
        raise NotImplementedError()

    def ___Pr_sgMC_nodes___(self, rp):
        raise NotImplementedError()

    def __call__(self, degree_plus):
        """

        Parameters
        ----------
        degree_plus : int
            the Gauss degree will locally be `cell.space.N + degree_plus`.

        Returns
        -------

        """
        assert degree_plus >= 0 and degree_plus % 1 == 0, f"degree_plus={degree_plus} is wrong."
        self._dp_ = degree_plus
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
            cell = self._mesh_[rp]
            N = cell.N

            if N in self._coo_:
                pass
            else:
                _1dNodes, _1dWeights = Quadrature(N + self._dp_, category='Gauss').quad
                _2dNodes = np.meshgrid(_1dNodes, _1dNodes, indexing='ij')
                # noinspection PyUnresolvedReferences
                _2dNodes_ravel = [_.ravel('F') for _ in _2dNodes]
                _2dWeights = np.tensordot(_1dWeights, _1dWeights, axes=0)
                _2dWeights_ravel = _2dWeights.ravel('F')

                self._coo_[N] = [_1dNodes  , _2dNodes  , _2dNodes_ravel], \
                                [_1dWeights, _2dWeights, _2dWeights_ravel]

            return self._coo_[N]

        #------------ segment -----------------------------------------------------------------
        elif rp.__class__.__name__ == 'mpRfT2_Segment':

            N = rp.N

            if N in self._coo_:
                pass
            else:
                _1dNodes, _1dWeights = Quadrature(N + self._dp_, category='Gauss').quad
                self._coo_[N] = _1dNodes, _1dWeights

            return self._coo_[N]

        #------------ else -----------------------------------------------------------------
        else:
            raise NotImplementedError(f'Not implemented Gauss coo_map for {rp}')
        #===================================================================================




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/coo_map/helpers/Gauss.py
    pass
