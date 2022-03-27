

from screws.freeze.base import FrozenOnly
import numpy as np


class _2dCSCG_So1F_Reconstruct(FrozenOnly):
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
        :param vectorized:
        :param value_only:
        :return:
        """
        f = self._f_
        mesh = self._f_.mesh

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

            assert INDICES == mesh.elements.indices, f"currently, vectorized computation only works" \
                                                          f"for full reconstruction."

            iJ = mesh.elements.coordinate_transformation.vectorized.inverse_Jacobian_matrix(*xietasigma)

            if len(INDICES) == 0:
                vx, vy = None, None
            else:
                numOfBasisComponents = f.num.basis_components
                Arr = f.cochain.array
                ArrX = Arr[:,:numOfBasisComponents[0]]
                ArrY = Arr[:,numOfBasisComponents[0]:]

                u = np.einsum('ij, ki -> kj', basis[0], ArrX, optimize='greedy')
                v = np.einsum('ij, ki -> kj', basis[1], ArrY, optimize='greedy')

                if mesh.elements.IS.homogeneous_according_to_types_wrt_metric:
                    if mesh.elements.IS.all_orthogonal:
                        vx = + np.einsum('kj, j -> kj', u, iJ[1][1], optimize='greedy')
                        vy = + np.einsum('kj, j -> kj', v, iJ[0][0], optimize='greedy')
                    else:
                        vx = + np.einsum('kj, j -> kj', u, iJ[1][1], optimize='greedy') \
                             - np.einsum('kj, j -> kj', v, iJ[0][1], optimize='greedy')
                        vy = - np.einsum('kj, j -> kj', u, iJ[1][0], optimize='greedy') \
                             + np.einsum('kj, j -> kj', v, iJ[0][0], optimize='greedy')

                else:
                    if mesh.elements.IS.all_orthogonal:
                        vx = + np.einsum('kj, kj -> kj', u, iJ[1][1], optimize='greedy')
                        vy = + np.einsum('kj, kj -> kj', v, iJ[0][0], optimize='greedy')
                    else:
                        vx = + np.einsum('kj, kj -> kj', u, iJ[1][1], optimize='greedy') \
                             - np.einsum('kj, kj -> kj', v, iJ[0][1], optimize='greedy')
                        vy = - np.einsum('kj, kj -> kj', u, iJ[1][0], optimize='greedy') \
                             + np.einsum('kj, kj -> kj', v, iJ[0][0], optimize='greedy')

            if ravel:
                pass
            else:
                raise NotImplementedError()

            if value_only:
                return [vx, vy]
            else:
                raise Exception()

        # ----- non-vectorized ------------------------------------------------
        else:
            if value_only:
                raise NotImplementedError()
            else:
                xyz = dict()
                value = dict()
                shape = [len(xi), len(eta)]
                iJ = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
                for i in INDICES:
                    element = mesh.elements[i]
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    u = np.einsum('ij, i -> j', basis[0], f.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
                    v = np.einsum('ij, i -> j', basis[1], f.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')
                    value[i] = [None, None]
                    typeWr2Metric = element.type_wrt_metric.mark
                    iJi = iJ[i]
                    if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                        value[i][0] = u*iJi[1][1]
                        value[i][1] = v*iJi[0][0]
                    else:
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
    