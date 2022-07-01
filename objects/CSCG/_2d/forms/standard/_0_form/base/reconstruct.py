# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly
import numpy as np


class _2dCSCG_S0F_Reconstruct(FrozenOnly):
    """"""
    def __init__(self, f):
        self._f_ = f
        self._freeze_self_()



    def __call__(self, xi, eta, ravel=False, i=None, vectorized=False, value_only=False):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param i:
        :param vectorized: {bool} Do we vectorize the computation?
        :param value_only: {bool} Do we only compute the values, no xyz?
        :return:
            - if vectorized: return a 3d array of the reconstructed values only.
            - if not vectorized: return dictionaries of coordinates and value.
        """
        f = self._f_
        mesh = self._f_.mesh

        xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta)

        #-------- parse indices --------------------------------------------------------
        if i is None:
            INDICES = mesh.elements.indices
        else:
            if vectorized:
                vectorized = False

            if isinstance(i, int):
                INDICES = [i,]
            else:
                raise Exception()
        # ----------- vectorized reconstruction ----------------------------------------
        if vectorized:

            assert INDICES == mesh.elements.indices, f"currently, vectorized computation only works" \
                                                          f"for full reconstruction."
            if len(INDICES) > 0:
                v = np.einsum('ij, ki -> kj', basis[0], f.cochain.array, optimize='greedy')

            else:
                v = None

            if ravel:
                pass
            else:
                raise NotImplementedError()

            if value_only:
                return (v,)
            else:
                raise Exception()

        # ----------- non-vectorized reconstruction ----------------------------------------
        else:

            if value_only:
                raise Exception()

            else:

                xyz = dict()
                value = dict()
                shape = [len(xi), len(eta)]
                for i in INDICES:
                    element = mesh.elements[i]
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    v = np.einsum('ij, i -> j', basis[0], f.cochain.local[i], optimize='greedy')
                    if ravel:
                        value[i] = [v,]
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                        value[i] = [v.reshape(shape, order='F'),]
                return xyz, value