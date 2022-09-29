# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 12:22 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np

class miUsTriangular_S2F_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, xi, eta, ravel=False, i=None, vectorized=False, value_only=False):
        """
        Reconstruct the standard 3-form.

        Given ``xi``, ``eta`` and ``sigma``, we reconstruct the 3-form on ``meshgrid(xi, eta, sigma)``
        in all elements.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param i: (`default`:``None``) Do the reconstruction for ``#i`` element. if it is ``None``,
            then do it for all elements.
        :type i: int, None
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :param bool ravel: (`default`:``False``) If we return 1d data?
        :param vectorized:
        :param value_only:
        :returns: A tuple of outputs

            1. (Dict[int, list]) -- :math:`x, y, z` coordinates.
            2. (Dict[int, list]) -- Reconstructed values.
        """
        f = self._sf_
        mesh = self._sf_.mesh

        xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta)

        #--- parse indices --------------------------------------------------
        if i is None: # default, in all local mesh-elements.
            INDICES = mesh.elements.indices
        else:
            if vectorized: vectorized = False

            if isinstance(i ,int):
                INDICES = [i, ]
            else:
                raise NotImplementedError()
        #---- vectorized -----------------------------------------------
        if vectorized:

            raise NotImplementedError()

            # assert INDICES == mesh.elements.indices, f"currently, vectorized computation only works" \
            #                                               f"for full reconstruction."
            #
            # det_iJ = mesh.elements.coordinate_transformation.vectorized.inverse_Jacobian(*xietasigma)
            #
            # if len(INDICES) > 0:
            #     if mesh.elements.IS.homogeneous_according_to_types_wrt_metric:
            #         v = np.einsum('ij, ki, j -> kj', basis[0], f.cochain.array, det_iJ, optimize='greedy')
            #     else:
            #         v = np.einsum('ij, ki, kj -> kj', basis[0], f.cochain.array, det_iJ, optimize='greedy')
            #
            # else:
            #     v = None
            #
            # if ravel:
            #     pass
            # else:
            #     raise NotImplementedError()
            #
            # if value_only:
            #     return v,
            # else:
            #     raise Exception()

        #----- non-vectorized ------------------------------------------------
        else:
            if value_only:
                raise NotImplementedError()
            else:
                xyz = dict()
                value = dict()
                shape = [len(xi), len(eta)]
                for i in INDICES:
                    element = mesh.elements[i]
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                    basis_det_iJ = basis[0] * det_iJ

                    v = np.einsum('ij, i -> j', basis_det_iJ, f.cochain.local[i], optimize='greedy')
                    if ravel:
                        value[i] = [v,]
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                        value[i] = [v.reshape(shape, order='F'),]
                return xyz, value


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
