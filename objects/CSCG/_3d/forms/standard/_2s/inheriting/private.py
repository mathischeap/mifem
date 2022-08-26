# -*- coding: utf-8 -*-


from scipy.sparse import csc_matrix, bmat
import numpy as np



# noinspection PyUnresolvedReferences
class _3dCSCG_S2F_Private:

    def ___PRIVATE_make_reconstruction_matrix_on_grid___(self, xi, et, sg, element_range=None):
        """Make the reconstruction matrices for all mesh elements. These matrices are stored in
        a dict whose keys are the numbers of mesh elements and values are the local reconstruction
        matrices.

        Let `RM` be the reconstruction matrix (or the tuple of three matrices).
        If we want to do the local reconstruction, we do

            RM[i] @ f.cochain.local[i]

        and we will get the reconstructions of the form `f` on `meshgrid(xi, eta, sigma)` in mesh-element
        #i. And if `f` is a scalar form, we get a 1d array. And if `f` is a vector form, we get a
        tuple of three 1d arrays (its three components along x, y, z directions.)

        :param xi: 1d array
        :param et: 1d array
        :param sg: 1d array
        :param element_range:
            We are going to construct matrices for these mesh elements. It can be one of
                1) None: for all local elements
                2) 'mesh boundary': those local elements are attached to mesh boundary.
        :return:
        """
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, et, sg)
        iJ = self.mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
        b0, b1, b2 = basis
        b0 = b0.T
        b1 = b1.T
        b2 = b2.T
        OO01 = 0 * b1
        OO02 = 0 * b2
        OO10 = 0 * b0
        OO12 = 0 * b2
        OO20 = 0 * b0
        OO21 = 0 * b1
        CACHE = dict()
        #------ take care of element range -----------------------------------------------
        if element_range is None:
            INDICES = self.mesh.elements
        elif element_range == 'mesh boundary':
            INDICES = self.mesh.boundaries.involved_elements
        else:
            raise Exception(f"element_range = {element_range} is wrong!")

        #-------- generating the reconstruction matrices for included mesh elements -------------
        RM = dict()
        for i in INDICES:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in CACHE:
                    RM[i] = CACHE[typeWr2Metric]
                else:
                    iJi = iJ[i]
                    _0u = iJi[1][1] * iJi[2][2] - iJi[1][2] * iJi[2][1]
                    _1v = iJi[2][2] * iJi[0][0] - iJi[2][0] * iJi[0][2]
                    _2w = iJi[0][0] * iJi[1][1] - iJi[0][1] * iJi[1][0]
                    rm00 = np.einsum('ji, j -> ji', b0, _0u, optimize='greedy')
                    rm11 = np.einsum('ji, j -> ji', b1, _1v, optimize='greedy')
                    rm22 = np.einsum('ji, j -> ji', b2, _2w, optimize='greedy')

                    if typeWr2Metric[:4] == 'Orth':
                        RM_i_ = ( np.hstack((rm00, OO01, OO02)),
                                  np.hstack((OO10, rm11, OO12)),
                                  np.hstack((OO20, OO21, rm22)) )
                    else:
                        _0v = iJi[2][1] * iJi[0][2] - iJi[2][2] * iJi[0][1]
                        _0w = iJi[0][1] * iJi[1][2] - iJi[0][2] * iJi[1][1]
                        _1u = iJi[1][2] * iJi[2][0] - iJi[1][0] * iJi[2][2]
                        _1w = iJi[0][2] * iJi[1][0] - iJi[0][0] * iJi[1][2]
                        _2u = iJi[1][0] * iJi[2][1] - iJi[1][1] * iJi[2][0]
                        _2v = iJi[2][0] * iJi[0][1] - iJi[2][1] * iJi[0][0]
                        rm01 = np.einsum('ji, j -> ji', b1, _0v, optimize='greedy')
                        rm02 = np.einsum('ji, j -> ji', b2, _0w, optimize='greedy')
                        rm10 = np.einsum('ji, j -> ji', b0, _1u, optimize='greedy')
                        rm12 = np.einsum('ji, j -> ji', b2, _1w, optimize='greedy')
                        rm20 = np.einsum('ji, j -> ji', b0, _2u, optimize='greedy')
                        rm21 = np.einsum('ji, j -> ji', b1, _2v, optimize='greedy')
                        RM_i_ = (np.hstack((rm00, rm01, rm02)),
                                 np.hstack((rm10, rm11, rm12)),
                                 np.hstack((rm20, rm21, rm22)))

                    CACHE[typeWr2Metric] = RM_i_
                    RM[i] = RM_i_
            else:
                iJi = iJ[i]
                _0u = iJi[1][1] * iJi[2][2] - iJi[1][2] * iJi[2][1]
                _0v = iJi[2][1] * iJi[0][2] - iJi[2][2] * iJi[0][1]
                _0w = iJi[0][1] * iJi[1][2] - iJi[0][2] * iJi[1][1]
                _1u = iJi[1][2] * iJi[2][0] - iJi[1][0] * iJi[2][2]
                _1v = iJi[2][2] * iJi[0][0] - iJi[2][0] * iJi[0][2]
                _1w = iJi[0][2] * iJi[1][0] - iJi[0][0] * iJi[1][2]
                _2u = iJi[1][0] * iJi[2][1] - iJi[1][1] * iJi[2][0]
                _2v = iJi[2][0] * iJi[0][1] - iJi[2][1] * iJi[0][0]
                _2w = iJi[0][0] * iJi[1][1] - iJi[0][1] * iJi[1][0]

                rm00 = np.einsum('ji, j -> ji', b0, _0u, optimize='greedy')
                rm01 = np.einsum('ji, j -> ji', b1, _0v, optimize='greedy')
                rm02 = np.einsum('ji, j -> ji', b2, _0w, optimize='greedy')

                rm10 = np.einsum('ji, j -> ji', b0, _1u, optimize='greedy')
                rm11 = np.einsum('ji, j -> ji', b1, _1v, optimize='greedy')
                rm12 = np.einsum('ji, j -> ji', b2, _1w, optimize='greedy')

                rm20 = np.einsum('ji, j -> ji', b0, _2u, optimize='greedy')
                rm21 = np.einsum('ji, j -> ji', b1, _2v, optimize='greedy')
                rm22 = np.einsum('ji, j -> ji', b2, _2w, optimize='greedy')

                RM[i] = (np.hstack((rm00, rm01, rm02)),
                         np.hstack((rm10, rm11, rm12)),
                         np.hstack((rm20, rm21, rm22)))

        return RM


    def ___PRIVATE_operator_inner___(self, other, i, xietasigma, quad_weights, bfSelf, bfOther):
        """
        We compute the inner product between ``self`` and ``other`` in element ``i``.

        Note that here we only return a local matrix.
        """
        element = self.mesh.elements[i]
        mark = element.type_wrt_metric.mark
        J = element.coordinate_transformation.Jacobian_matrix(*xietasigma)
        sqrtg = element.coordinate_transformation.Jacobian(*xietasigma, J=J)
        iJ = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=J)
        g = element.coordinate_transformation.inverse_metric_matrix(*xietasigma, iJ=iJ)
        del J, iJ

        if isinstance(mark, str) and mark[:4] == 'Orth':
            M00 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[1][1]*g[2][2], bfOther[0], bfSelf[0])
            M11 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[2][2]*g[0][0], bfOther[1], bfSelf[1])
            M22 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[0][0]*g[1][1], bfOther[2], bfSelf[2])
            M01 = None
            M02 = None
            M12 = None
            M10 = None
            M20 = None
            M21 = None

        else:
            M00 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[1][1]*g[2][2]-g[1][2]*g[2][1], bfOther[0], bfSelf[0])
            M11 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[2][2]*g[0][0]-g[2][0]*g[0][2], bfOther[1], bfSelf[1])
            M22 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[0][0]*g[1][1]-g[0][1]*g[1][0], bfOther[2], bfSelf[2])
            g12_20_g10_22 = g[1][2] * g[2][0] - g[1][0] * g[2][2]
            g10_21_g11_20 = g[1][0] * g[2][1] - g[1][1] * g[2][0]
            g20_01_g21_00 = g[2][0] * g[0][1] - g[2][1] * g[0][0]
            M01 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g12_20_g10_22, bfOther[0], bfSelf[1])
            M02 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g10_21_g11_20, bfOther[0], bfSelf[2])
            M12 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g20_01_g21_00, bfOther[1], bfSelf[2])
            if other is self:
                M10 = M01.T
                M20 = M02.T
                M21 = M12.T
            else:
                M10 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g12_20_g10_22, bfOther[1], bfSelf[0])
                M20 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g10_21_g11_20, bfOther[2], bfSelf[0])
                M21 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g20_01_g21_00, bfOther[2], bfSelf[1])

        Mi = bmat([(M00, M01, M02),
                   (M10, M11, M12),
                   (M20, M21, M22)], format='csc')
        return Mi

    @staticmethod
    def ___PRIVATE_inner_H1___(quad_weights, sqrt_g, g, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights * sqrt_g * g, bfO, bfS, optimize='optimal')
        return csc_matrix(M)
