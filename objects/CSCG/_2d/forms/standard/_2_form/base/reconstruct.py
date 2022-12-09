# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')
from objects.CSCG._2d.forms.standard.base.reconstruct import _2dCSCG_SF_ReconstructBase
from objects.CSCG._2d.discreteFields.scalar.main import _2dCSCG_DF_Scalar
import numpy as np


class _2dCSCG_S2F_Reconstruct(_2dCSCG_SF_ReconstructBase):
    """"""
    def __init__(self, f):
        super(_2dCSCG_S2F_Reconstruct, self).__init__(f)
        self._freeze_self_()


    def __call__(self, xi, eta, ravel=False, element_range=None, vectorized=False, value_only=False):
        """
        Reconstruct the 2d standard 2-form.

        Given ``xi``, ``eta`` and ``sigma``, we reconstruct the 3-form on ``meshgrid(xi, eta, sigma)``
        in all elements.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param element_range: (`default`:``None``) Do the reconstruction for ``#i`` element. if it is ``None``,
            then do it for all elements.
        :type element_range: int, None
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :param bool ravel: (`default`:``False``) If we return 1d data?
        :param vectorized:
        :param value_only:
        :returns: A tuple of outputs

            1. (Dict[int, list]) -- :math:`x, y, z` coordinates.
            2. (Dict[int, list]) -- Reconstructed values.
        """
        f = self._f_
        mesh = self._f_.mesh

        xietasigma, basis = f.do.evaluate_basis_at_meshgrid(xi, eta)
        #--- parse indices --------------------------------------------------
        if element_range is None: # default, in all local mesh-elements.
            INDICES = mesh.elements.indices
        else:
            if vectorized: vectorized = False

            if isinstance(element_range, int):
                INDICES = [element_range,]
            else:
                raise NotImplementedError()

        #---- vectorized -----------------------------------------------
        if vectorized:

            assert INDICES == mesh.elements.indices, f"currently, vectorized computation only works" \
                                                          f"for full reconstruction."

            det_iJ = mesh.elements.coordinate_transformation.vectorized.inverse_Jacobian(*xietasigma)

            if len(INDICES) > 0:
                if mesh.elements.whether.homogeneous_according_to_types_wrt_metric:
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

        #----- non-vectorized ------------------------------------------------
        else:
            if value_only:
                raise NotImplementedError()
            else:
                xyz = dict()
                value = dict()
                shape = [len(xi), len(eta)]
                iJC = dict()
                for i in INDICES:
                    element = mesh.elements[i]
                    typeWr2Metric = element.type_wrt_metric.mark
                    xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
                    if typeWr2Metric in iJC:
                        basis_det_iJ = iJC[typeWr2Metric]
                    else:
                        det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                        basis_det_iJ = basis[0] * det_iJ
                        if isinstance(typeWr2Metric, str):
                            iJC[typeWr2Metric] = basis_det_iJ

                    v = np.einsum('ij, i -> j', basis_det_iJ, f.cochain.local[i], optimize='greedy')
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
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                v = np.einsum('ij, i -> j', basis[0] * det_iJ, f.cochain.local[e], optimize='greedy')

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