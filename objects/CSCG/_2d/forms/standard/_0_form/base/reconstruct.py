# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')
from objects.CSCG._2d.forms.standard.base.reconstruct import _2dCSCG_SF_ReconstructBase
from objects.CSCG._2d.discreteFields.scalar.main import _2dCSCG_DF_Scalar
import numpy as np


class _2dCSCG_S0F_Reconstruct(_2dCSCG_SF_ReconstructBase):
    """"""
    def __init__(self, f):
        super(_2dCSCG_S0F_Reconstruct, self).__init__(f)
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
                return v, # do not remove comma
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
                    # noinspection PyUnresolvedReferences
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    v = np.einsum('ij, i -> j', basis[0], f.cochain.local[i], optimize='greedy')
                    if ravel:
                        value[i] = [v,]
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                        value[i] = [v.reshape(shape, order='F'),]

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

        #-------- reconstructing -----------------------------------------------------------------1
        xy = dict()
        value = dict()
        EMPTY_DATA = np.empty((0,0))

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
                #___DIFF for different forms____ reconstruction in local mesh element #e________diff
                xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, et)
                xy[e] = element.coordinate_transformation.mapping(*xietasigma)
                v = np.einsum('ij, i -> j', basis[0], f.cochain.local[e], optimize='greedy')

                # noinspection PyUnresolvedReferences
                xy[e] = [xy[e][j].reshape(shape, order='F') for j in range(2)]
                # noinspection PyUnresolvedReferences
                value[e] = [v.reshape(shape, order='F'),] #=============================diff

        #-------- prime-region-wise stack coordinates and values ----------------------------------1
        XY, VAL, element_global_numbering = self.___PRIVATE_distribute_XY_and_VAL___(mesh, xy, value)

        XY = self.___PRIVATE_prime_region_wise_stack___(mesh, XY, 2, grid, element_global_numbering)
        VAL = self.___PRIVATE_prime_region_wise_stack___(mesh, VAL, 1, grid, element_global_numbering)

        return _2dCSCG_DF_Scalar(mesh, XY, VAL,
                                 name=self._f_.standard_properties.name,
                                 structured=True, grid=grid)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/_0_form/base/reconstruct.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.1)([30,30])
    # mesh = MeshGenerator('chp1',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)
    ES = ExactSolutionSelector(mesh)('sL:sincos1')
    f0 = FC('0-f-o', is_hybrid=True)
    f0.TW.func.do.set_func_body_as(ES, 'potential')

    f0.TW.current_time = 0
    f0.TW.do.push_all_to_instant()
    f0.discretize()

    x = np.linspace(-1,1,10)

    ds = f0.reconstruct.discrete_scalar([x, x])
    ds.visualize()