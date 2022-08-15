# -*- coding: utf-8 -*-

from screws.freeze.base import FrozenOnly

import numpy as np


class _3dCSCG_SF2_reconstruct(FrozenOnly):
    """"""
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel=False, i=None, regions=None, vectorized=False, value_only=False):
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