# -*- coding: utf-8 -*-


import sys
if './' not in sys.path: sys.path.append('./')


from objects.CSCG._3d.forms.standard.base.reconstruct import _3dCSCG_SF_Reconstruct
import numpy as np
from objects.CSCG._3d.discreteDields.scalar.main import _3dCSCG_DF_Scalar






class _3dCSCG_SF3_Reconstruct(_3dCSCG_SF_Reconstruct):
    """"""
    def __init__(self, sf):
        super(_3dCSCG_SF3_Reconstruct, self).__init__(sf)
        self._freeze_self_()


    def __call__(self, xi, eta, sigma, ravel=False, i=None, regions=None, vectorized=False, value_only=False):
        """
        Reconstruct the standard 3-form.

        Given ``xi``, ``eta`` and ``sigma``, we reconstruct the 3-form on ``meshgrid(xi, eta, sigma)``
        in all elements.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param sigma: A 1d iterable object of floats between -1 and 1.
        :param i: (`default`:``None``) Do the reconstruction for ``#i`` element. if it is ``None``,
            then do it for all elements.
        :type i: int, None
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :type sigma: list, tuple, numpy.ndarray
        :param bool ravel: (`default`:``False``) If we return 1d data?
        :param regions: Higher priority than input ``i``.
        :param vectorized:
        :param value_only:
        :returns: A tuple of outputs

            1. (Dict[int, list]) -- :math:`x, y, z` coordinates.
            2. (Dict[int, list]) -- Reconstructed values.
        """
        f = self._sf_
        mesh = f.mesh

        xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta, sigma)

        #----------- parse INDICES -----------------------------------------------------------------
        if regions is None:
            if i is None:
                INDICES = mesh.elements.indices
            elif isinstance(i, int):
                INDICES = [i, ]
            else:
                raise Exception(f"i={i} is wrong.")
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
                ri = mesh.do.FIND_region_name_of_element(i)
                if ri in regions:
                    INDICES.append(i)

        #---- vectorized ---------------------------------------------------------------------------
        if vectorized:

            assert INDICES == mesh.elements.indices, f"currently, vectorized computation only works" \
                                                          f"for full reconstruction."

            det_iJ = mesh.elements.coordinate_transformation.vectorized.inverse_Jacobian(*xietasigma)

            if len(INDICES) > 0:

                if mesh.elements.IS.homogeneous_according_to_types_wrt_metric:

                    v = np.einsum('ij, ki, j -> kj', basis[0], f.cochain.array, det_iJ, optimize='greedy')

                else:

                    v = np.einsum('ij, ki, kj -> kj', basis[0], f.cochain.array, det_iJ, optimize='greedy')

            else:
                v = None

            if ravel:
                pass
            else:
                raise NotImplementedError()

            if value_only:
                return v,
            else:
                raise Exception()

        #-------- non-vectorized -------------------------------------------------------------------
        else:
            shape = [len(xi), len(eta), len(sigma)]
            value = dict()
            if value_only:
                raise Exception()
            else:
                xyz = dict()
                biJC = dict()
                for i in INDICES:
                    element = mesh.elements[i]
                    typeWr2Metric = element.type_wrt_metric.mark
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    if typeWr2Metric in biJC:
                        basis_det_iJ = biJC[typeWr2Metric]
                    else:
                        det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                        basis_det_iJ = basis[0] * det_iJ
                        if isinstance(typeWr2Metric, str):
                            biJC[typeWr2Metric] = basis_det_iJ
                    v = np.einsum('ij, i -> j', basis_det_iJ, f.cochain.local[i], optimize='greedy')
                    if ravel:
                        value[i] = [v,]
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                        value[i] = [v.reshape(shape, order='F'),]
                return xyz, value
        #===========================================================================================


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
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)

                xyz[e] = element.coordinate_transformation.mapping(*xietasigma)
                v = np.einsum('ij, i -> j', basis[0] * det_iJ, f.cochain.local[e], optimize='optimal')
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
    # mpiexec -n 5 python objects/CSCG/_3d/forms/standard/_3s/reconstruct.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('bridge_arch_cracked')([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t

    scalar = FC('scalar', u)

    f = FC('3-f', is_hybrid=False)
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