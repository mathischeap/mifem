
# -*- coding: utf-8 -*-

import numpy as np
from objects.CSCG._3d.forms.standard.base.reconstruct import _3dCSCG_SF_Reconstruct


class _3dCSCG_SF0_reconstruct(_3dCSCG_SF_Reconstruct):
    """"""
    def __init__(self, sf):
        super(_3dCSCG_SF0_reconstruct, self).__init__(sf)
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

        xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta, sigma)

        # ---- parse INDICES ----------------------------------------------------------------------
        if regions is None:
            if i is None:
                INDICES = mesh.elements.indices
            elif isinstance(i, int):
                INDICES = [i, ]
            else:
                raise Exception()
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

        # ----------- vectorized reconstruction --------------------------------------------------
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
                return (v,)
            else:
                raise Exception()

        #------- non-vectorized -------------------------------------------------------------------
        else:

            if value_only:
                raise NotImplementedError()
            else:

                xyz = dict()
                value = dict()
                shape = [len(xi), len(eta), len(sigma)]
                for i in INDICES:
                    element = mesh.elements[i]
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    v = np.einsum('ij, i -> j', basis[0], f.cochain.local[i], optimize='optimal')
                    if ravel:
                        value[i] = [v,]
                    else:
                        # noinspection PyUnresolvedReferences
                        xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                        value[i] = [v.reshape(shape, order='F'),]
                return xyz, value