# -*- coding: utf-8 -*-

import sys
if './' not in sys.path:
    sys.path.append('./')

import numpy as np
from objects.CSCG._2d.forms.standard.base.reconstruct import _2dCSCG_SF_ReconstructBase
from objects.CSCG._2d.discreteFields.vector.main import _2dCSCG_DF_Vector


class _2dCSCG_So1F_Reconstruct(_2dCSCG_SF_ReconstructBase):
    """"""
    def __init__(self, f):
        super(_2dCSCG_So1F_Reconstruct, self).__init__(f)
        self._freeze_self_()

    def __call__(self, xi, eta, ravel=False, element_range=None, vectorized=False, value_only=False):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param element_range:
        :param vectorized:
        :param value_only:
        :return:
        """
        f = self._f_
        mesh = self._f_.mesh

        xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta)

        # --- parse indices --------------------------------------------------
        if element_range is None:  # default, in all local mesh-elements.
            INDICES = mesh.elements.indices
        else:
            if vectorized:
                vectorized = False

            if isinstance(element_range, int):
                INDICES = [element_range, ]
            else:
                raise NotImplementedError()

        # ---- vectorized -----------------------------------------------
        if vectorized:

            assert INDICES == mesh.elements.indices, f"currently, vectorized computation only works" \
                                                          f"for full reconstruction."

            iJ = mesh.elements.coordinate_transformation.vectorized.inverse_Jacobian_matrix(*xietasigma)

            if len(INDICES) == 0:
                vx, vy = None, None
            else:
                numOfBasisComponents = f.num.basis_components
                Arr = f.cochain.array
                ArrX = Arr[:, :numOfBasisComponents[0]]
                ArrY = Arr[:, numOfBasisComponents[0]:]

                u = np.einsum('ij, ki -> kj', basis[0], ArrX, optimize='greedy')
                v = np.einsum('ij, ki -> kj', basis[1], ArrY, optimize='greedy')

                if mesh.elements.whether.homogeneous_according_to_types_wrt_metric:
                    if mesh.elements.whether.all_orthogonal:
                        vx = + np.einsum('kj, j -> kj', u, iJ[1][1], optimize='greedy')
                        vy = + np.einsum('kj, j -> kj', v, iJ[0][0], optimize='greedy')

                    else:
                        vx = + np.einsum(
                            'kj, j -> kj', u, iJ[1][1], optimize='greedy'
                        ) - np.einsum(
                            'kj, j -> kj', v, iJ[0][1], optimize='greedy'
                        )

                        vy = - np.einsum(
                            'kj, j -> kj', u, iJ[1][0], optimize='greedy'
                        ) + np.einsum(
                            'kj, j -> kj', v, iJ[0][0], optimize='greedy'
                        )

                else:
                    if mesh.elements.whether.all_orthogonal:
                        vx = + np.einsum('kj, kj -> kj', u, iJ[1][1], optimize='greedy')
                        vy = + np.einsum('kj, kj -> kj', v, iJ[0][0], optimize='greedy')
                    else:
                        vx = + np.einsum(
                            'kj, kj -> kj', u, iJ[1][1], optimize='greedy'
                        ) - np.einsum(
                            'kj, kj -> kj', v, iJ[0][1], optimize='greedy'
                        )

                        vy = - np.einsum(
                            'kj, kj -> kj', u, iJ[1][0], optimize='greedy'
                        ) + np.einsum(
                            'kj, kj -> kj', v, iJ[0][0], optimize='greedy'
                        )

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
                xyz: dict = dict()
                value: dict = dict()
                shape: list = [len(xi), len(eta)]
                iJ = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
                for i in INDICES:
                    element = mesh.elements[i]
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    u = np.einsum('ij, i -> j',
                                  basis[0], f.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
                    v = np.einsum('ij, i -> j',
                                  basis[1], f.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')
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

    def discrete_vector(self, grid):
        """We reconstruct this outer S1F as a 2d CSCG discrete vector in the `regions`.

        Parameters
        ----------
        grid: dict, list, tuple
            For example:
                rst = [r, s, ], these r, s will be used for all regions.
                rst = {'R:R1': [r1, s1], 'R:R2': [r2, s2], ....}, in each
                    region, we use its own r, s

            r or s can be in [-1,1] or [0,1].

        Returns
        -------

        """
        f = self._f_
        mesh = f.mesh
        Xi_Eta_Sigma_D, grid = self.___PRIVATE_distribute_region_wise_meshgrid___(mesh, grid)

        # -------- reconstructing -----------------------------------------------------------------1
        xy = dict()
        value = dict()
        EMPTY_DATA = np.empty((0, 0))

        for e in mesh.elements:
            rn, ij = mesh.do.find.region_name_and_local_indices_of_element(e)
            X, E = Xi_Eta_Sigma_D[rn]
            i, j = ij
            # these xi, et, sg are for reconstruction in local mesh element [e]
            xi, et = X[i], E[j]
            element = mesh.elements[e]
            shape = [len(xi), len(et)]

            if shape[0] < 1 or shape[1] < 1:
                xy[e] = [EMPTY_DATA, EMPTY_DATA]
                value[e] = [EMPTY_DATA, EMPTY_DATA]
            else:
                # ___DIFF for different forms____ reconstruction in local mesh element #e________diff
                xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, et)
                iJi = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
                xy[e] = element.coordinate_transformation.mapping(*xietasigma)
                u = np.einsum('ij, i -> j', basis[0], f.cochain.___PRIVATE_local_on_axis___('x', e), optimize='greedy')
                v = np.einsum('ij, i -> j', basis[1], f.cochain.___PRIVATE_local_on_axis___('y', e), optimize='greedy')

                value[e] = [None, None]
                typeWr2Metric = element.type_wrt_metric.mark

                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    value[e][0] = u*iJi[1][1]
                    value[e][1] = v*iJi[0][0]
                else:
                    value[e][0] = + u*iJi[1][1] - v*iJi[0][1]
                    value[e][1] = - u*iJi[1][0] + v*iJi[0][0]

                # noinspection PyUnresolvedReferences
                xy[e] = [xy[e][j].reshape(shape, order='F') for j in range(2)]
                # noinspection PyUnresolvedReferences
                value[e] = [value[e][j].reshape(shape, order='F') for j in range(2)]  # =========diff

        # -------- prime-region-wise stack coordinates and values ----------------------------------1
        XY, VAL, element_global_numbering = self.___PRIVATE_distribute_XY_and_VAL___(mesh, xy, value)

        XY = self.___PRIVATE_prime_region_wise_stack___(mesh, XY, 2, grid, element_global_numbering)
        VAL = self.___PRIVATE_prime_region_wise_stack___(mesh, VAL, 2, grid, element_global_numbering)

        return _2dCSCG_DF_Vector(mesh, XY, VAL,
                                 name=self._f_.standard_properties.name,
                                 structured=True, grid=grid)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/_1_form/outer/reconstruct.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    # mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    mesh = MeshGenerator('chp2')([10, 10])
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)

    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    u = FC('1-f-o', is_hybrid=True)

    u.CF = ES.velocity
    u.CF.current_time = 0
    u.discretize()

    r = np.linspace(-1, 1, 5)
    s = np.linspace(-1, 1, 6)

    dv = u.reconstruct.discrete_vector([r, s])
