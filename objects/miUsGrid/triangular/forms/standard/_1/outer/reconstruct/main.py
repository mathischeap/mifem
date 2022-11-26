# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 7:08 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
import numpy as np

class miUsTriangular_oS1F_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()


    def __call__(self, xi, eta, ravel=False, i=None, vectorized=False, value_only=False):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param i: The elements in which we want to reconstruct.
        :param vectorized: {bool} Do we vectorize the computation?
        :param value_only: {bool} Do we only compute the values, no xyz?
        :return:
            - if vectorized: return a 3d array of the reconstructed values only.
            - if not vectorized: return dictionaries of coordinates and value.
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

        # ----- non-vectorized ------------------------------------------------
        else:
            value = dict()
            shape = [len(xi), len(eta)]
            IJI = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)

            for i in INDICES:

                u = np.einsum('ij, i -> j', basis[i][0], f.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
                v = np.einsum('ij, i -> j', basis[i][1], f.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')

                value[i] = [None, None]
                iJi = IJI[i]

                value[i][0] = + u*iJi[1][1] - v*iJi[0][1]
                value[i][1] = - u*iJi[1][0] + v*iJi[0][0]

                if ravel:
                    pass
                else:
                    # noinspection PyUnresolvedReferences
                    value[i] = [value[i][j].reshape(shape, order='F') for j in range(2)]

            if value_only:
                return value
            else:
                xyz = dict()

                for i in INDICES:
                    element = mesh.elements[i]
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)

                    if ravel:
                        pass
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]

                return xyz, value


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
