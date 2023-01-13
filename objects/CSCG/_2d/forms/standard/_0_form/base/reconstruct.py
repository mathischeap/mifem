# -*- coding: utf-8 -*-
import sys
if './' not in sys.path:
    sys.path.append('./')
from objects.CSCG._2d.forms.standard.base.reconstruct import _2dCSCG_SF_ReconstructBase
from objects.CSCG._2d.discreteFields.scalar.main import _2dCSCG_DF_Scalar
import numpy as np


class _2dCSCG_S0F_Reconstruct(_2dCSCG_SF_ReconstructBase):
    """"""
    def __init__(self, f):
        super(_2dCSCG_S0F_Reconstruct, self).__init__(f)
        self._freeze_self_()

    def __call__(self, xi, eta, ravel=False, element_range=None, vectorized=False, value_only=False):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param element_range:
        :param vectorized: {bool} Do we vectorize the computation?
        :param value_only: {bool} Do we only compute the values, no xyz?
        :return:
            - if vectorized: return a 3d array of the reconstructed values only.
            - if not vectorized: return dictionaries of coordinates and value.
        """
        f = self._f_
        mesh = self._f_.mesh

        xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta)

        # -------- parse indices --------------------------------------------------------
        if element_range is None:
            INDICES = mesh.elements.indices

        elif isinstance(element_range, str):
            if vectorized:
                vectorized = False
            else:
                pass

            bns = mesh.boundaries.names
            if element_range in bns:

                INDICES = mesh.boundaries.range_of_element_edges[element_range]

            else:
                raise NotImplementedError(f"can not parse element_range = {element_range}.")

        else:
            if vectorized:
                vectorized = False
            else:
                pass

            if isinstance(element_range, int):
                INDICES = [element_range, ]
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
                return v,  # do not remove comma
            else:
                raise Exception()

        # ----------- non-vectorized reconstruction ----------------------------------------
        else:

            if value_only:
                raise Exception()

            else:

                xyz: dict = dict()
                value: dict = dict()
                shape: list = [len(xi), len(eta)]

                for ind in INDICES:
                    if isinstance(ind, int):
                        i = ind
                        position = None
                    elif isinstance(ind, str):  # we get a str indicator
                        if ind[:-1].isnumeric() and ind[-1] in 'UDLR':  # the str indicator indicates a element edge.
                            i = int(ind[:-1])
                            position = ind[-1]
                        else:
                            raise NotImplementedError()
                    else:
                        raise NotImplementedError()

                    element = mesh.elements[i]

                    if position is None:  # when position is None, we reconstruct all over the element.
                        xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                        v = np.einsum('ij, i -> j', basis[0], f.cochain.local[i], optimize='greedy')
                        if ravel:
                            value[i] = [v, ]
                        else:
                            # noinspection PyUnresolvedReferences
                            xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                            value[i] = [v.reshape(shape, order='F'), ]

                    elif position in ['U', 'D', 'L', 'R']:  # we will reconstruct only on the edge indicated by position

                        # the algorithm below is not optimal, but it is fine for 2d.

                        if position == 'U':
                            assert xi[0] == -1, f"to construct on U edge, xi must be [-1, ...]"
                        elif position == 'D':
                            assert xi[-1] == 1, f"to construct on D edge, xi must be [..., 1]"
                        elif position == 'L':
                            assert eta[0] == 1, f"to construct on L edge, eta must be [-1, ...]"
                        elif position == 'R':
                            assert eta[-1] == 1, f"to construct on R edge, xi must be [..., 1]"
                        else:
                            raise Exception()

                        coo = element.coordinate_transformation.mapping(*xietasigma)
                        v = np.einsum('ij, i -> j', basis[0], f.cochain.local[i], optimize='greedy')

                        x, y = [coo[j].reshape(shape, order='F') for j in range(2)]
                        v = v.reshape(shape, order='F')

                        if position == 'U':
                            x, y, v = x[0, :], y[0, :], v[0, :]
                        elif position == 'D':
                            x, y, v = x[-1, :], y[-1, :], v[-1, :]
                        elif position == 'L':
                            x, y, v = x[:, 0], y[:, 0], v[:, 0]
                        elif position == 'R':
                            x, y, v = x[:, -1], y[:, -1], v[:, -1]
                        else:
                            raise Exception()

                        if i not in xyz:
                            xyz[i] = dict()
                            value[i] = dict()
                        else:
                            pass

                        xyz[i][position] = (x, y)
                        value[i][position] = (v, )

                    else:
                        raise NotImplementedError(f"position={position} not implemented.")

                return xyz, value

    def discrete_scalar(self, grid):
        """We reconstruct this outer S0F as a 2d CSCG discrete scalar in the `regions`.

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
                # ___DIFF for different forms____ reconstruction in local mesh element # e ________diff
                xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, et)
                xy[e] = element.coordinate_transformation.mapping(*xietasigma)
                v = np.einsum('ij, i -> j', basis[0], f.cochain.local[e], optimize='greedy')

                # noinspection PyUnresolvedReferences
                xy[e] = [xy[e][j].reshape(shape, order='F') for j in range(2)]
                # noinspection PyUnresolvedReferences
                value[e] = [v.reshape(shape, order='F'), ]  # =============================diff

        # -------- prime-region-wise stack coordinates and values ----------------------------------1
        XY, VAL, element_global_numbering = self.___PRIVATE_distribute_XY_and_VAL___(mesh, xy, value)

        XY = self.___PRIVATE_prime_region_wise_stack___(mesh, XY, 2, grid, element_global_numbering)
        VAL = self.___PRIVATE_prime_region_wise_stack___(mesh, VAL, 1, grid, element_global_numbering)

        return _2dCSCG_DF_Scalar(mesh, XY, VAL,
                                 name=self._f_.standard_properties.name,
                                 structured=True, grid=grid)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/_0_form/base/reconstruct.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.1)([5, 5])

    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)
    ES = ExactSolutionSelector(mesh)('sL:sincos1')
    f0 = FC('0-f-o', hybrid=True)
    f0.CF = ES.potential
    f0.CF.current_time = 0
    f0.discretize()

    x = np.linspace(-1, 1, 100)

    ds = f0.reconstruct.discrete_scalar([x, x])
    ds.visualize()
