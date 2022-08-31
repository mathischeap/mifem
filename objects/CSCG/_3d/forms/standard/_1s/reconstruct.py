# -*- coding: utf-8 -*-


import numpy as np
from objects.CSCG._3d.forms.standard.base.reconstruct import _3dCSCG_SF_Reconstruct


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