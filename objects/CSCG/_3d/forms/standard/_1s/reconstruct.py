# -*- coding: utf-8 -*-


import sys
if './' not in sys.path: sys.path.append('./')


from root.config.main import np
from objects.CSCG._3d.forms.standard.base.reconstruct import _3dCSCG_SF_Reconstruct
from objects.CSCG._3d.discreteDields.vector.main import _3dCSCG_DF_Vector


class _3dCSCG_SF1_reconstruct(_3dCSCG_SF_Reconstruct):
    """"""
    def __init__(self, sf):
        super(_3dCSCG_SF1_reconstruct, self).__init__(sf)
        self._freeze_self_()





    def __call__(self, xi, eta, sigma, ravel=False, i=None, regions=None, vectorized=False, value_only=False):
        """

        :param xi:
        :param eta:
        :param sigma:
        :param ravel:
        :param i:
            In which elements we do the reconstruction?
        :param regions: Higher priority than input i.
        :param vectorized:
        :param value_only:
        :return:
        """
        f = self._sf_
        mesh = self._sf_.mesh

        #------ parse INDICES -------------------------------------------------------------
        if regions is None:
            if i is None:
                INDICES = mesh.elements.indices
            elif isinstance(i, int):
                INDICES = [i, ]
            else:
                raise NotImplementedError()
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
        #------------ vectorized -----------------------------------------------------------------
        if vectorized:
            raise NotImplementedError()

        # ------- non-vectorized -----------------------------------------------------------------
        else:
            if value_only:
                raise NotImplementedError()
            else:
                xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta, sigma)
                xyz = dict()
                value = dict()
                shape = [len(xi), len(eta), len(sigma)]

                iJ = mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
                for i in INDICES:
                    element = mesh.elements[i]
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    u = np.einsum('ij, i -> j', basis[0], f.cochain.___PRIVATE_local_on_axis___('x', i), optimize='optimal')
                    v = np.einsum('ij, i -> j', basis[1], f.cochain.___PRIVATE_local_on_axis___('y', i), optimize='optimal')
                    w = np.einsum('ij, i -> j', basis[2], f.cochain.___PRIVATE_local_on_axis___('z', i), optimize='optimal')
                    value[i] = [None, None, None]
                    typeWr2Metric = element.type_wrt_metric.mark
                    iJi = iJ[i]
                    if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                        value[i][0] = u*iJi[0][0]
                        value[i][1] = v*iJi[1][1]
                        value[i][2] = w*iJi[2][2]
                    else:
                        for j in range(3):
                            value[i][j] = u*iJi[0][j] + v*iJi[1][j] + w*iJi[2][j]
                    if ravel:
                        pass
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                        # noinspection PyUnresolvedReferences
                        value[i] = [value[i][j].reshape(shape, order='F') for j in range(3)]
                return xyz, value


    def discrete_vector(self, grid):
        """We reconstruct this S2F as a 3d CSCG discrete vector in the `regions`.

        Parameters
        ----------
        grid: dict, list, tuple
            For example:
                rst = [r, s, t], these r, s, t will be used for all regions.
                rst = {'R:R1': [r1, s1, t1], 'R:R2': [r2, s2, t2], ....}, in each
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
                iJi = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
                xyz[e] = element.coordinate_transformation.mapping(*xietasigma)
                u = np.einsum('ij, i -> j', basis[0], f.cochain.___PRIVATE_local_on_axis___('x', e), optimize='greedy')
                v = np.einsum('ij, i -> j', basis[1], f.cochain.___PRIVATE_local_on_axis___('y', e), optimize='greedy')
                w = np.einsum('ij, i -> j', basis[2], f.cochain.___PRIVATE_local_on_axis___('z', e), optimize='greedy')

                value[e] = [None, None, None]
                typeWr2Metric = element.type_wrt_metric.mark

                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                        value[e][0] = u*iJi[0][0]
                        value[e][1] = v*iJi[1][1]
                        value[e][2] = w*iJi[2][2]

                else:
                    for j in range(3):
                        value[e][j] = u*iJi[0][j] + v*iJi[1][j] + w*iJi[2][j]

                # noinspection PyUnresolvedReferences
                xyz[e] = [xyz[e][j].reshape(shape, order='F') for j in range(3)]
                # noinspection PyUnresolvedReferences
                value[e] = [value[e][j].reshape(shape, order='F') for j in range(3)] #=========diff

        #-------- prime-region-wise stack coordinates and values ----------------------------------1
        XYZ, VAL, element_global_numbering = self.___PRIVATE_distribute_XYZ_and_VAL___(mesh, xyz, value)
        XYZ = self.___PRIVATE_prime_region_wise_stack___(mesh, XYZ, 3, grid, element_global_numbering)
        VAL = self.___PRIVATE_prime_region_wise_stack___(mesh, VAL, 3, grid, element_global_numbering)

        return _3dCSCG_DF_Vector(mesh, XYZ, VAL, name=self._sf_.standard_properties.name,
                                 structured=True, grid=grid)





if __name__ == '__main__':
    # mpiexec -n 5 python objects/CSCG/_3d/forms/standard/_1s/reconstruct.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('bridge_arch_cracked')([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    velocity = FC('vector', (u,v,w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)

    f1 = FC('1-f', is_hybrid=False)
    f1.TW.func.do.set_func_body_as(velocity)
    f1.TW.current_time = 0
    f1.TW.___DO_push_all_to_instant___()
    f1.discretize()
    # f2.visualize(x=0.25)

    r = [0,0.25, 0.5, 0.75, 1]
    s = [0,0.125,0.25, 0.5, 0.75, 1]
    t = [0,0.25, 0.5,0.6, 0.75, 1]

    dv = f1.reconstruct.discrete_vector([r, s, t])

    print(dv.coordinates)