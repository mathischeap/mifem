# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 7:08 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
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

            # assert INDICES == mesh.elements.indices, f"currently, vectorized computation only works" \
            #                                               f"for full reconstruction."
            #
            # iJ = mesh.elements.coordinate_transformation.vectorized.inverse_Jacobian_matrix(*xietasigma)
            #
            # if len(INDICES) == 0:
            #     vx, vy = None, None
            # else:
            #     numOfBasisComponents = f.num.basis_components
            #     Arr = f.cochain.array
            #     ArrX = Arr[:,:numOfBasisComponents[0]]
            #     ArrY = Arr[:,numOfBasisComponents[0]:]
            #
            #     u = np.einsum('ij, ki -> kj', basis[0], ArrX, optimize='greedy')
            #     v = np.einsum('ij, ki -> kj', basis[1], ArrY, optimize='greedy')
            #
            #     if mesh.elements.IS.homogeneous_according_to_types_wrt_metric:
            #         if mesh.elements.IS.all_orthogonal:
            #             vx = + np.einsum('kj, j -> kj', u, iJ[1][1], optimize='greedy')
            #             vy = + np.einsum('kj, j -> kj', v, iJ[0][0], optimize='greedy')
            #         else:
            #             vx = + np.einsum('kj, j -> kj', u, iJ[1][1], optimize='greedy') \
            #                  - np.einsum('kj, j -> kj', v, iJ[0][1], optimize='greedy')
            #             vy = - np.einsum('kj, j -> kj', u, iJ[1][0], optimize='greedy') \
            #                  + np.einsum('kj, j -> kj', v, iJ[0][0], optimize='greedy')
            #
            #     else:
            #         if mesh.elements.IS.all_orthogonal:
            #             vx = + np.einsum('kj, kj -> kj', u, iJ[1][1], optimize='greedy')
            #             vy = + np.einsum('kj, kj -> kj', v, iJ[0][0], optimize='greedy')
            #         else:
            #             vx = + np.einsum('kj, kj -> kj', u, iJ[1][1], optimize='greedy') \
            #                  - np.einsum('kj, kj -> kj', v, iJ[0][1], optimize='greedy')
            #             vy = - np.einsum('kj, kj -> kj', u, iJ[1][0], optimize='greedy') \
            #                  + np.einsum('kj, kj -> kj', v, iJ[0][0], optimize='greedy')
            #
            # if ravel:
            #     pass
            # else:
            #     raise NotImplementedError()
            #
            # if value_only:
            #     return [vx, vy]
            # else:
            #     raise Exception()

        # ----- non-vectorized ------------------------------------------------
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

                    u = np.einsum('ij, i -> j', basis[0], f.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
                    v = np.einsum('ij, i -> j', basis[1], f.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')


                    value[i] = [None, None]
                    iJi = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)

                    value[i][0] = + u*iJi[1][1] - v*iJi[0][1]
                    value[i][1] = - u*iJi[1][0] + v*iJi[0][0]

                    if ravel:
                        pass
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                        # noinspection PyUnresolvedReferences
                        value[i] = [value[i][j].reshape(shape, order='F') for j in range(2)]

                return xyz, value


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
