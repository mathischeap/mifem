# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:55 PM
"""
from components.freeze.main import FrozenOnly
import numpy as np
from scipy.sparse import csr_matrix, bmat

class _3dCSCG_2LocalTrace_Reconstruct(FrozenOnly):
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
        element_range : mesh element numbers

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

        xietasigma, pb = self._ltf_.do.evaluate_basis_at_meshgrid(xi, eta, sigma)
        ii, jj, kk = np.size(xi), np.size(eta), np.size(sigma)
        xyz : dict = dict()
        v : dict = dict()

        Tmap = mesh.trace.elements.map
        for i in indices:
            if i in mesh.elements:

                tes = Tmap[i]
                xyz_i : list = list()
                v_i : list = list()

                for te_number, side in zip(tes, 'NSWEBF'):
                    te = mesh.trace.elements[te_number]

                    if i in self._ltf_.cochain.local:

                        side_prime_cochain = \
                            self._ltf_.cochain.___components_of_cochain_on_element_side___(
                                i, side
                            )

                    else:

                        if i in self._ltf_.cochain.local_ESW:
                            esw_i = self._ltf_.cochain.local_ESW[i]
                            if side in esw_i:
                                side_prime_cochain = esw_i[side]
                            else:
                                side_prime_cochain = None
                        else:
                            side_prime_cochain = None

                    if side_prime_cochain is not None:

                        xyz_i.append(
                            te.coordinate_transformation.mapping(
                                *xietasigma[side], from_element=i, side=side
                            )
                        )
                        g = te.coordinate_transformation.metric(*xietasigma[side])

                        if side in 'NS':
                            vi = np.einsum('i, j, ij -> j', side_prime_cochain, 1 / np.sqrt(g),
                                           pb[side][0], optimize='greedy')
                        elif side in 'WE':
                            vi = np.einsum('i, j, ij -> j', side_prime_cochain, 1 / np.sqrt(g),
                                           pb[side][0], optimize='greedy')
                        elif side in 'BF':
                            vi = np.einsum('i, j, ij -> j', side_prime_cochain, 1 / np.sqrt(g),
                                           pb[side][0], optimize='greedy')
                        else:
                            raise Exception()

                        v_i.append(vi)

                    else:
                        xyz_i.append(None)
                        v_i.append(None)

                if all([_ is None for _ in v_i]):
                    continue
                else:
                    pass

                XYZ = dict()
                V = dict()

                if ravel:

                    for j, side in enumerate('NSWEBF'):
                        XYZ[side] = xyz_i[j]
                        V[side] = [v_i[j],]

                else:
                    XYZ = dict()
                    V = dict()
                    for j, side in enumerate('NSWEBF'):
                        xyz_ = xyz_i[j]
                        v_ = v_i[j]

                        if xyz_ is None:
                            pass
                        else:
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
        cacheDict : dict = dict()

        Tmap = mesh.trace.elements.map
        for i in indices:
            if i in mesh.elements:
                element = mesh.elements[i]
                mark = element.type_wrt_metric.mark

                if isinstance(mark, str) and mark in cacheDict:

                    RD[i], RD_sides[i] = cacheDict[mark]

                else:

                    tes = Tmap[i]
                    M = list()

                    for te_number, side in zip(tes, 'NSWEBF'):
                        te = mesh.trace.elements[te_number]
                        g = te.coordinate_transformation.metric(*xietasigma[side])

                        if side in 'NS':
                            ms = np.einsum('j, ij -> ji', 1 / np.sqrt(g),
                                           pb[side][0], optimize='greedy')
                        elif side in 'WE':
                            ms = np.einsum('j, ij -> ji', 1 / np.sqrt(g),
                                           pb[side][0], optimize='greedy')
                        elif side in 'BF':
                            ms = np.einsum('j, ij -> ji', 1 / np.sqrt(g),
                                           pb[side][0], optimize='greedy')
                        else:
                            raise Exception()

                        M.append(ms)

                    RDi = bmat(
                        [
                            [csr_matrix(M[0]), None, None, None, None, None],
                            [None, csr_matrix(M[1]), None, None, None, None],
                            [None, None, csr_matrix(M[2]), None, None, None],
                            [None, None, None, csr_matrix(M[3]), None, None],
                            [None, None, None, None, csr_matrix(M[4]), None],
                            [None, None, None, None, None, csr_matrix(M[5])],
                        ],
                        format='csr'
                    )

                    RDi_sides = {
                        'N': M[0],
                        'S': M[1],
                        'W': M[2],
                        'E': M[3],
                        'B': M[4],
                        'F': M[5],
                    }

                    if isinstance(mark, str):
                        cacheDict[mark] = (RDi, RDi_sides)
                    else:
                        pass

                    RD[i] = RDi
                    RD_sides[i] = RDi_sides

        return RD, RD_sides