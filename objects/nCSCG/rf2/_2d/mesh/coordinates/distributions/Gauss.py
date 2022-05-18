# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/17/2022 12:18 AM
"""
import sys

import numpy as np

if './' not in sys.path: sys.path.append('./')
# import numpy as np
from screws.quadrature import Quadrature

from objects.nCSCG.rf2._2d.mesh.coordinates.distributions.base import _2nCSCG_MRF2_CooDistributionBase


class Gauss(_2nCSCG_MRF2_CooDistributionBase):
    """"""

    def __init__(self, mesh):
        """"""
        super(Gauss, self).__init__(mesh)
        self._distribution_ = 'Gauss'
        self._dp_ = None
        self._coo_ = dict()
        self._freeze_self_()

    def __call__(self, degree_plus):
        """

        Parameters
        ----------
        degree_plus : int
            the Gauss degree will locally be `cell.space.N + degree_plus`.
        ndim :
            We return `ndim`-dimensional coordinates.

        Returns
        -------

        """
        assert self.signature == self._mesh_.signature, f"signature dis-match"
        assert degree_plus >= 0 and degree_plus % 1 == 0, f"degree_plus={degree_plus} is wrong."
        self._dp_ = degree_plus
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
        cell = self._mesh_(indices)
        N = cell.space.N

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




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
