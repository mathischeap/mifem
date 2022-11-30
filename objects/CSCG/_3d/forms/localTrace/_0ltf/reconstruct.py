# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:55 PM
"""
from components.freeze.main import FrozenOnly

import numpy as np
from scipy.sparse import csr_matrix, bmat

class _3dCSCG_0LocalTrace_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel=False, element_range=None):
        """

        Parameters
        ----------
        xi
        eta
        sigma
        ravel
        element_range : mesh element number

        Returns
        -------

        """
        mesh = self._ltf_.mesh

        if element_range is None:
            indices = mesh.elements._elements_.keys()
        elif isinstance(element_range, int):
            indices = [element_range,]
        else:
            indices = element_range

        xietasigma, pb = self._ltf_.do.evaluate_basis_at_meshgrid(xi, eta,sigma)
        ii, jj, kk = np.size(xi), np.size(eta), np.size(sigma)
        xyz = dict()
        v = dict()

        Tmap = mesh.trace.elements.map
        for i in indices:
            if i in mesh.elements:

                tes = Tmap[i]
                xyz_i = list()
                v_i = list()
                for te_number, side in zip(tes, 'NSWEBF'):
                    te = mesh.trace.elements[te_number]

                    xyz_i.append(
                        te.coordinate_transformation.mapping(
                            *xietasigma[side], from_element=i, side=side
                        )
                    )

                    side_prime_cochain = self._ltf_.cochain.___components_of_cochain_on_element_side___(
                        i, side
                    )

                    v_i.append(
                        np.einsum(
                            'i, ij -> j',
                            side_prime_cochain,
                            pb[side][0],
                            optimize='greedy'
                        )
                    )

                if ravel:
                    X, Y, Z = list(), list(), list()

                    for xyz_ in xyz_i:
                        x, y, z = xyz_
                        X.append(x)
                        Y.append(y)
                        Z.append(z)

                    xyz[i] = [
                        np.concatenate(X),
                        np.concatenate(Y),
                        np.concatenate(Z)
                    ]

                    v[i] = [np.concatenate(v_i),]

                else:
                    XYZ = dict()
                    V = dict()
                    for j, side in enumerate('NSWEBF'):
                        xyz_ = xyz_i[j]
                        v_ = v_i[j]
                        if side in 'NS':
                            xyz_ = [xyz_[_].reshape((jj, kk), order='F') for _ in range(3)]
                            v_ = [v_.reshape((jj, kk), order='F'),]
                        elif side in 'WE':
                            xyz_ = [xyz_[_].reshape((ii, kk), order='F') for _ in range(3)]
                            v_ = [v_.reshape((ii, kk), order='F'), ]
                        elif side in 'BF':
                            xyz_ = [xyz_[_].reshape((ii, jj), order='F') for _ in range(3)]
                            v_ = [v_.reshape((ii, jj), order='F'), ]
                        else:
                            raise Exception

                        XYZ[side] = xyz_
                        V[side] = v_

                    xyz[i] = XYZ
                    v[i] = V

        return xyz, v

    def ___PrLT_make_reconstruction_matrix_on_grid___(self, xi, et, sg, element_range=None):
        """

        Parameters
        ----------
        xi
        et
        sg
        element_range

        Returns
        -------

        """
        mesh = self._ltf_.mesh

        if element_range is None:
            indices = mesh.elements._elements_.keys()
        elif isinstance(element_range, int):
            indices = [element_range,]
        else:
            indices = element_range

        xietasigma, pb = self._ltf_.do.evaluate_basis_at_meshgrid(xi, et, sg)
        RD = dict()
        RD_sides = dict()

        M = list()
        for side in 'NSWEBF':
            M.append(csr_matrix(pb[side][0].T))

        rdi = bmat(
            [
                [M[0], None, None, None, None, None],
                [None, M[1], None, None, None, None],
                [None, None, M[2], None, None, None],
                [None, None, None, M[3], None, None],
                [None, None, None, None, M[4], None],
                [None, None, None, None, None, M[5]],
            ],
            format='csr'
        )

        rdi_sides = {
            'N': M[0],
            'S': M[1],
            'W': M[2],
            'E': M[3],
            'B': M[4],
            'F': M[5],
        }

        for i in indices:
            if i in mesh.elements:
                RD[i] = rdi
                RD_sides[i] = rdi_sides

        return RD, RD_sides