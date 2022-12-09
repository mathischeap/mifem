# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')

import numpy as np
from objects.CSCG._3d.forms.standard.base.reconstruct import _3dCSCG_SF_Reconstruct
from objects.CSCG._3d.discreteDields.scalar.main import _3dCSCG_DF_Scalar


class _3dCSCG_SF0_reconstruct(_3dCSCG_SF_Reconstruct):
    """"""
    def __init__(self, sf):
        super(_3dCSCG_SF0_reconstruct, self).__init__(sf)
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel=False, element_range=None, regions=None, vectorized=False, value_only=False):
        """

        :param xi:
        :param eta:
        :param sigma:
        :param ravel:
        :param element_range:
            In which elements we do the reconstruction?
        :param regions: Higher priority than input i.
        :param vectorized:
        :param value_only:
        :return:
        """
        f = self._sf_
        mesh = self._sf_.mesh

        xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta, sigma)

        # ---- parse INDICES -----------------------------------------------------------
        if regions is None:
            if element_range is None:
                INDICES = mesh.elements.indices
            elif isinstance(element_range, int):
                INDICES = [element_range,]
            else:
                INDICES = element_range
        else:
            if regions == 'all':
                regions = mesh.domain.regions
            elif isinstance(regions, str):
                regions = [regions,]
            else:
                pass
            assert isinstance(regions, (list, tuple)), f"regions={regions} is wrong."
            assert len(set(regions)) == len(regions), f"regions={regions} has repeated regions."
            for i, r in enumerate(regions):
                assert r in mesh.domain.regions, f"regions[{i}]={r} is wrong."

            INDICES = list()
            for i in mesh.elements.indices:
                ri = mesh.do.find.region_name_of_element(i)
                if ri in regions:
                    INDICES.append(i)

        # ----------- vectorized reconstruction ----------------------------------------
        if vectorized:

            assert INDICES == mesh.elements.indices, \
                f"currently, vectorized computation only works " \
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

        #------- non-vectorized --------------------------------------------------------
        else:

            if value_only:
                raise NotImplementedError()
            else:

                xyz = dict()
                value = dict()
                shape = [len(xi), len(eta), len(sigma)]
                for i in INDICES:
                    element = mesh.elements[i]
                    # noinspection PyUnresolvedReferences
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    v = np.einsum('ij, i -> j', basis[0], f.cochain.local[i], optimize='optimal')
                    if ravel:
                        value[i] = [v,]
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                        value[i] = [v.reshape(shape, order='F'),]
                return xyz, value

    def discrete_scalar(self, grid):
        """We reconstruct this S2F as a 3d CSCG discrete vector in the `regions`.

        Parameters
        ----------
        grid: dict, list, tuple
            For example:
                grid = [r, s, t], these r, s, t will be used for all regions.
                grid = {'R:R1': [r1, s1, t1], 'R:R2': [r2, s2, t2], ....}, in each
                    region, we use its own r, s, t

        Returns
        -------

        """
        f = self._sf_
        mesh = f.mesh
        Xi_Eta_Sigma_D, grid = self.___PRIVATE_distribute_region_wise_meshgrid___(mesh, grid)
        EMPTY_DATA = np.empty((0,0,0))

        #-------- reconstructing -----------------------------------------------------------------1
        xyz = dict()
        value = dict()

        for e in mesh.elements:
            rn, ijk = mesh.do.find.region_name_and_local_indices_of_element(e)
            X, E, S = Xi_Eta_Sigma_D[rn]
            i, j, k = ijk
            # these xi, et, sg are for reconstruction in local mesh element [e]
            xi, et, sg = X[i], E[j], S[k]
            element = mesh.elements[e]
            shape = [len(xi), len(et), len(sg)]

            if shape[0] < 1 or shape[1] < 1 or shape[2] < 1:
                xyz[e] = [EMPTY_DATA, EMPTY_DATA, EMPTY_DATA]
                value[e] = [EMPTY_DATA, EMPTY_DATA, EMPTY_DATA]
            else:
                #___DIFF for different forms____ reconstruction in local mesh element #e________diff
                xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, et, sg)
                xyz[e] = element.coordinate_transformation.mapping(*xietasigma)
                v = np.einsum('ij, i -> j', basis[0], f.cochain.local[e], optimize='optimal')
                # noinspection PyUnresolvedReferences
                xyz[e] = [xyz[e][j].reshape(shape, order='F') for j in range(3)]
                # noinspection PyUnresolvedReferences
                value[e] = [v.reshape(shape, order='F'),] #=========diff

        #-------- prime-region-wise stack coordinates and values ----------------------------------1
        XYZ, VAL, element_global_numbering = self.___PRIVATE_distribute_XYZ_and_VAL___(mesh, xyz, value)
        XYZ = self.___PRIVATE_prime_region_wise_stack___(mesh, XYZ, 3, grid, element_global_numbering)
        VAL = self.___PRIVATE_prime_region_wise_stack___(mesh, VAL, 1, grid, element_global_numbering)

        return _3dCSCG_DF_Scalar(mesh, XYZ, VAL, name=self._sf_.standard_properties.name,
                                 structured=True, grid=grid)

if __name__ == '__main__':
    # mpiexec -n 5 python objects/CSCG/_3d/forms/standard/_0s/reconstruct.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('bridge_arch_cracked')([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t

    scalar = FC('scalar', u)

    f = FC('0-f', is_hybrid=False)
    f.TW.func.do.set_func_body_as(scalar)
    f.TW.current_time = 0
    f.TW.___DO_push_all_to_instant___()
    f.discretize()
    # f2.visualize(x=0.25)

    r = [0,0.25, 0.5, 0.75, 1]
    s = [0,0.125,0.25, 0.5, 0.75, 1]
    t = [0,0.25, 0.5,0.6, 0.75, 1]

    dv = f.reconstruct.discrete_scalar([r, s, t])

    print(dv.coordinates)