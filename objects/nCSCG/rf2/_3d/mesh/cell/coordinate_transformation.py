# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from screws.freeze.base import FrozenOnly

class _3nCSCG_CellCT(FrozenOnly):
    """"""
    def __init__(self, cell):
        self._cell_ = cell
        self._me_ct_ = self._cell_.mesh.cscg.elements[self._cell_.indices[0]].coordinate_transformation
        self._freeze_self_()


    @property
    def origin_and_delta(self):
        """consider the level-0 cell (cscg mesh-element) as a [-1,1]^3 reference domain.
        """
        return self._cell_.mesh.do.find.origin_and_delta(self._cell_.indices)

    def ___PRIVATE_plot_data___(self, zoom=1, density=None):
        """

        Parameters
        ----------
        zoom
        density

        Returns
        -------
        data : ndarray
            of shape (12, 3, density) which stands for 12 lines, xyz coordinates and density.

        """
        if density is None:
            mark = self._cell_.type_wrt_metric.mark
            if isinstance(mark, str) and mark[:4] == 'Orth':
                density = 2
            else:
                density = 5
        else:
            pass

        assert 0 < zoom <= 1, f"zoom={zoom} is wrong!"

        O = np.ones(density) * zoom
        M = - np.ones(density) * zoom
        S = np.linspace(-1 * zoom, 1 * zoom, density)

        data = np.array([
        np.array(self.mapping(S, M, M)),
        np.array(self.mapping(S, M, O)),
        np.array(self.mapping(S, O, O)),
        np.array(self.mapping(S, O, M)),

        np.array(self.mapping(M, S, M)),
        np.array(self.mapping(M, S, O)),
        np.array(self.mapping(O, S, O)),
        np.array(self.mapping(O, S, M)),

        np.array(self.mapping(M, M, S)),
        np.array(self.mapping(M, O, S)),
        np.array(self.mapping(O, O, S)),
        np.array(self.mapping(O, M, S))
        ])

        return data

    def mapping(self, xi, et, sg):
        """

        Parameters
        ----------
        xi :
            Any dimension, in [-1,1]. xi, et, sg same shape
        et :
            Any dimension, in [-1,1]. xi, et, sg same shape
        sg :
            Any dimension, in [-1,1]. xi, et, sg same shape


        Returns
        -------

        """
        o, d = self.origin_and_delta
        xi = o[0] + (xi + 1) * d / 2
        et = o[1] + (et + 1) * d / 2
        sg = o[2] + (sg + 1) * d / 2
        return self._me_ct_.mapping(xi, et, sg)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_3d/mesh/cell/coordinate_transformation.py
    from objects.nCSCG.rf2._3d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None)
    i = 8
    if i in mesh.cscg.elements:

        cell = mesh((i, 3, 3), dynamic=True)

        ct = cell.coordinate_transformation

        xi = np.linspace(-1,1,5)
        eta = np.linspace(-1,1,8)
        sigma = np.linspace(-1,1,8)
        xi, eta, sigma = np.meshgrid(xi, eta, sigma, indexing='ij')

        xyz = ct.mapping(xi, eta, sigma)

        print(xyz)