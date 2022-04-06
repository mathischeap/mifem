



from screws.freeze.main import FrozenOnly
import numpy as np






class _3dCSCG_SF3_Reconstruct(FrozenOnly):
    def __init__(self, _3sf):
        self._sf_ = _3sf
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
                return (v,)
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