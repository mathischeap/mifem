# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')


from root.config.main import np
from objects.CSCG._3d.forms.standard.base.reconstruct import _3dCSCG_SF_Reconstruct
from objects.CSCG._3d.discrete_fields.vector.main import _3dCSCG_DF_Vector

class _3dCSCG_SF2_reconstruct(_3dCSCG_SF_Reconstruct):
    """"""
    def __init__(self, sf):
        super(_3dCSCG_SF2_reconstruct, self).__init__(sf)
        self._freeze_self_()

    def __call__(self, xi, eta, sigma,
                 ravel=False, i=None, regions=None,
                 vectorized=False, value_only=False):
        """
        Do the reconstruction.

        :param xi:
        :param eta:
        :param sigma:
        :param ravel:
        :param i: The element we want to reconstruct in. When it is ``None``, we do the reconstruction for
            all elements and store the results in one coordinate dictionary and one value dictionary.
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
                JC = dict() # local cache
                for i in INDICES:
                    element = mesh.elements[i]
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)

                    u = np.einsum('ij, i -> j', basis[0], f.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
                    v = np.einsum('ij, i -> j', basis[1], f.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')
                    w = np.einsum('ij, i -> j', basis[2], f.cochain.___PRIVATE_local_on_axis___('z', i), optimize='greedy')

                    value[i] = [None, None, None]
                    typeWr2Metric = element.type_wrt_metric.mark

                    if typeWr2Metric in JC:
                        _0u, _0v, _0w, _1u, _1v, _1w, _2u, _2v, _2w = JC[typeWr2Metric]

                    else:

                        iJi = iJ[i]

                        if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                            _0u = iJi[1][1] * iJi[2][2]
                            _0v = None
                            _0w = None
                            _1u = None
                            _1v = iJi[2][2] * iJi[0][0]
                            _1w = None
                            _2u = None
                            _2v = None
                            _2w = iJi[0][0] * iJi[1][1]
                        else:
                            _0u = iJi[1][1] * iJi[2][2] - iJi[1][2] * iJi[2][1]
                            _0v = iJi[2][1] * iJi[0][2] - iJi[2][2] * iJi[0][1]
                            _0w = iJi[0][1] * iJi[1][2] - iJi[0][2] * iJi[1][1]
                            _1u = iJi[1][2] * iJi[2][0] - iJi[1][0] * iJi[2][2]
                            _1v = iJi[2][2] * iJi[0][0] - iJi[2][0] * iJi[0][2]
                            _1w = iJi[0][2] * iJi[1][0] - iJi[0][0] * iJi[1][2]
                            _2u = iJi[1][0] * iJi[2][1] - iJi[1][1] * iJi[2][0]
                            _2v = iJi[2][0] * iJi[0][1] - iJi[2][1] * iJi[0][0]
                            _2w = iJi[0][0] * iJi[1][1] - iJi[0][1] * iJi[1][0]

                        if isinstance(typeWr2Metric, str):
                            JC[typeWr2Metric] = _0u, _0v, _0w, _1u, _1v, _1w, _2u, _2v, _2w

                    if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                        value[i][0] = u * _0u
                        value[i][1] = v * _1v
                        value[i][2] = w * _2w
                    else:
                        value[i][0] = u * _0u +  v * _0v +  w * _0w
                        value[i][1] = u * _1u +  v * _1v +  w * _1w
                        value[i][2] = u * _2u +  v * _2v +  w * _2w

                    if ravel:
                        pass
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                        # noinspection PyUnresolvedReferences
                        value[i] = [value[i][j].reshape(shape, order='F') for j in range(3)]
                return xyz, value

    def discrete_vector(self, rst):
        """We reconstruct this S2F as a 3d CSCG discrete vector in the `regions`.

        Parameters
        ----------
        rst: dict, list, tuple
            For example:
                rst = [r, s, t], these r, s, t will be used for all regions.
                rst = {'R:R1': [r1, s1, t1], 'R:R2': [r2, s2, t2], ....}, in each
                    region, we use its own r, s, t

        Returns
        -------

        """
        Xi_Eta_Sigma_D, rst = self.___PRIVATE_distribute_region_wise_meshgrid___(rst)
        f = self._sf_
        mesh = f.mesh
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
                    _0u = iJi[1][1] * iJi[2][2]
                    _1v = iJi[2][2] * iJi[0][0]
                    _2w = iJi[0][0] * iJi[1][1]

                    value[e][0] = u * _0u
                    value[e][1] = v * _1v
                    value[e][2] = w * _2w

                else:
                    _0u = iJi[1][1] * iJi[2][2] - iJi[1][2] * iJi[2][1]
                    _0v = iJi[2][1] * iJi[0][2] - iJi[2][2] * iJi[0][1]
                    _0w = iJi[0][1] * iJi[1][2] - iJi[0][2] * iJi[1][1]
                    _1u = iJi[1][2] * iJi[2][0] - iJi[1][0] * iJi[2][2]
                    _1v = iJi[2][2] * iJi[0][0] - iJi[2][0] * iJi[0][2]
                    _1w = iJi[0][2] * iJi[1][0] - iJi[0][0] * iJi[1][2]
                    _2u = iJi[1][0] * iJi[2][1] - iJi[1][1] * iJi[2][0]
                    _2v = iJi[2][0] * iJi[0][1] - iJi[2][1] * iJi[0][0]
                    _2w = iJi[0][0] * iJi[1][1] - iJi[0][1] * iJi[1][0]

                    value[e][0] = u * _0u + v * _0v + w * _0w
                    value[e][1] = u * _1u + v * _1v + w * _1w
                    value[e][2] = u * _2u + v * _2v + w * _2w

                # noinspection PyUnresolvedReferences
                xyz[e] = [xyz[e][j].reshape(shape, order='F') for j in range(3)]
                # noinspection PyUnresolvedReferences
                value[e] = [value[e][j].reshape(shape, order='F') for j in range(3)] #=========diff

        #-------- prime-region-wise stack coordinates and values ----------------------------------1
        XYZ, VAL, element_global_numbering = self.___PRIVATE_distribute_XYZ_and_VAL___(xyz, value)
        XYZ = self.___PRIVATE_prime_region_wise_stack___(XYZ, 3, rst, element_global_numbering)
        VAL = self.___PRIVATE_prime_region_wise_stack___(VAL, 3, rst, element_global_numbering)

        return _3dCSCG_DF_Vector(mesh, XYZ, VAL, name=self._sf_.standard_properties.name)




if __name__ == '__main__':
    # mpiexec -n 5 python objects/CSCG/_3d/forms/standard/_2s/reconstruct.py
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

    f2 = FC('2-f', is_hybrid=False)
    f2.TW.func.do.set_func_body_as(velocity)
    f2.TW.current_time = 0
    f2.TW.___DO_push_all_to_instant___()
    f2.discretize()
    # f2.visualize(x=0.25)

    r = [0,0.25, 0.5, 0.75, 1]
    s = [0,0.125,0.25, 0.5, 0.75, 1]
    t = [0,0.25, 0.5,0.6, 0.75, 1]

    dv = f2.reconstruct.discrete_vector([r, s, t])

    print(dv.regions)