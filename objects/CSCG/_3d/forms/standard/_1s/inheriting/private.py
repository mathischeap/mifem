# -*- coding: utf-8 -*-
from components.quadrature import Quadrature
import numpy as np
from scipy.sparse import csc_matrix, bmat, lil_matrix
from tools.linearAlgebra.dataStructures.globalMatrix.main import GlobalVector


# noinspection PyUnresolvedReferences
class _3dCSCG_S1F_Private:


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
        RM = dict()

        #------ take care of element range -----------------------------------------------
        if element_range is None:
            INDICES = self.mesh.elements
        elif element_range == 'mesh boundary':# all mesh-elements that are attached to any mesh boundary
            INDICES = self.mesh.boundaries.involved_elements
        else:
            raise Exception(f"element_range = {element_range} is wrong!")

        #-------- generating the reconstruction matrices for included mesh elements -------
        for i in INDICES:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in CACHE:
                    RM[i] = CACHE[typeWr2Metric]
                else:
                    iJi = iJ[i]
                    rm00 = np.einsum('ji, j -> ji', b0, iJi[0][0], optimize='greedy')
                    rm11 = np.einsum('ji, j -> ji', b1, iJi[1][1], optimize='greedy')
                    rm22 = np.einsum('ji, j -> ji', b2, iJi[2][2], optimize='greedy')

                    if typeWr2Metric[:4] == 'Orth':
                        RM_i_ = ( np.hstack((rm00, OO01, OO02)),
                                  np.hstack((OO10, rm11, OO12)),
                                  np.hstack((OO20, OO21, rm22)) )
                    else:
                        rm01 = np.einsum('ji, j -> ji', b1, iJi[1][0], optimize='greedy')
                        rm02 = np.einsum('ji, j -> ji', b2, iJi[2][0], optimize='greedy')
                        rm10 = np.einsum('ji, j -> ji', b0, iJi[0][1], optimize='greedy')
                        rm12 = np.einsum('ji, j -> ji', b2, iJi[2][1], optimize='greedy')
                        rm20 = np.einsum('ji, j -> ji', b0, iJi[0][2], optimize='greedy')
                        rm21 = np.einsum('ji, j -> ji', b1, iJi[1][2], optimize='greedy')
                        RM_i_ = ( np.hstack((rm00, rm01, rm02)),
                                  np.hstack((rm10, rm11, rm12)),
                                  np.hstack((rm20, rm21, rm22)))

                    CACHE[typeWr2Metric] = RM_i_
                    RM[i] = RM_i_
            else:
                iJi = iJ[i]

                rm00 = np.einsum('ji, j -> ji', b0, iJi[0][0], optimize='greedy')
                rm01 = np.einsum('ji, j -> ji', b1, iJi[1][0], optimize='greedy')
                rm02 = np.einsum('ji, j -> ji', b2, iJi[2][0], optimize='greedy')

                rm10 = np.einsum('ji, j -> ji', b0, iJi[0][1], optimize='greedy')
                rm11 = np.einsum('ji, j -> ji', b1, iJi[1][1], optimize='greedy')
                rm12 = np.einsum('ji, j -> ji', b2, iJi[2][1], optimize='greedy')

                rm20 = np.einsum('ji, j -> ji', b0, iJi[0][2], optimize='greedy')
                rm21 = np.einsum('ji, j -> ji', b1, iJi[1][2], optimize='greedy')
                rm22 = np.einsum('ji, j -> ji', b2, iJi[2][2], optimize='greedy')

                RM[i] = ( np.hstack((rm00, rm01, rm02)),
                          np.hstack((rm10, rm11, rm12)),
                          np.hstack((rm20, rm21, rm22)))

        return RM




    def ___PRIVATE_operator_inner___(self, other, i, xietasigma, quad_weights, bfSelf, bfOther):
        """Note that here we only return a local matrix."""
        element = self.mesh.elements[i]
        mark = element.type_wrt_metric.mark
        J = element.coordinate_transformation.Jacobian_matrix(*xietasigma)
        sqrtg = element.coordinate_transformation.Jacobian(*xietasigma, J=J)
        iJ = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=J)
        g = element.coordinate_transformation.inverse_metric_matrix(*xietasigma, iJ=iJ)
        del J, iJ
        M00 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[0][0], bfOther[0], bfSelf[0])
        M11 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[1][1], bfOther[1], bfSelf[1])
        M22 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[2][2], bfOther[2], bfSelf[2])

        if isinstance(mark, str) and mark[:4] == 'Orth':
            M01 = None
            M02 = None
            M12 = None
            M10 = None
            M20 = None
            M21 = None
        else:
            M01 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[0][1], bfOther[0], bfSelf[1])
            M02 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[0][2], bfOther[0], bfSelf[2])
            M12 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[1][2], bfOther[1], bfSelf[2])
            if other is self:
                M10 = M01.T
                M20 = M02.T
                M21 = M12.T
            else:
                M10 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[1][0], bfOther[1], bfSelf[0])
                M20 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[2][0], bfOther[2], bfSelf[0])
                M21 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[2][1], bfOther[2], bfSelf[1])

        Mi = bmat([(M00, M01, M02),
                   (M10, M11, M12),
                   (M20, M21, M22)], format='csc')

        return Mi

    @staticmethod
    def ___PRIVATE_inner_H1___(quad_weights, sqrt_g_g, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights*sqrt_g_g, bfO, bfS, optimize='optimal')
        return csc_matrix(M)


    # def ___PRIVATE_basis_boundary_integral_over_region_side___(self,
    #     vector, region_name, side_name, quad_degree=None):
    #     """
    #     TODO: to be completed.
    #
    #     :param vector:
    #     :param region_name:
    #     :param side_name:
    #     :param quad_degree:
    #     :return:
    #     """
    #     p = [self.dqp[i] + 3 for i in range(self.ndim)] if quad_degree is None else quad_degree
    #     quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
    #
    #     if side_name == 'F':
    #         quad_weights_2d = np.kron(quad_weights[1], quad_weights[0])
    #
    #
    #     v0, v1 = vector
    #     if isinstance(v0, (int, float)):
    #         pass
    #     else:
    #         raise NotImplementedError()
    #     if isinstance(v1, (int, float)):
    #         pass
    #     else:
    #         raise NotImplementedError()
    #
    #
    #     bfn_xi = self.space.basises[0].node_basis(x=quad_nodes[0])
    #     bfn_et = self.space.basises[1].node_basis(x=quad_nodes[1])
    #
    #
    #     bfe_xi = self.space.basises[0].edge_basis(x=quad_nodes[0])
    #     bfe_et = self.space.basises[1].edge_basis(x=quad_nodes[1])
    #
    #     if side_name == 'F':
    #         B = (np.kron(bfn_et, bfe_xi), np.kron(bfe_et, bfn_xi))
    #
    #     B0, B1 = B
    #
    #     elements = self.mesh.do.FIND_element_attach_to_region_side(region_name, side_name)
    #     elements = elements.ravel('F')
    #     ELEMENTS = list()
    #     for i in elements:
    #         if i in self.mesh.elements:
    #             ELEMENTS.append(i)
    #
    #
    #     GV = lil_matrix((1, self.GLOBAL_num_dofs))
    #
    #     for i in ELEMENTS:
    #         dofs = self.numbering.do.find.dofs_on_element_side(i, side_name)
    #         mark = self.mesh.elements[i].type_wrt_metric.mark
    #         spacing = self.mesh.elements[i].spacing
    #
    #         if mark[:4] == 'Orth':
    #             Lx, Ly, Lz = (spacing[0][1] - spacing[0][0]) * 2, \
    #                          (spacing[1][1] - spacing[1][0]) * 2, \
    #                          (spacing[2][1] - spacing[2][0]) * 2,
    #
    #             # print(spacing, flush=True)
    #             Jx = Lx / 2
    #             Jy = Ly / 2
    #             Jz = Lz / 2
    #
    #             if side_name in 'NS':
    #                 J0, J1 = Jy, Jz
    #             elif side_name in 'WE':
    #                 J0, J1 = Jx, Jz
    #             elif side_name in 'BF':
    #                 J0, J1 = Jx, Jy
    #             else:
    #                 raise Exception()
    #
    #             if isinstance(v0, (int, float)):
    #                 # noinspection PyUnboundLocalVariable
    #                 cewti1 = np.einsum('ki, i -> k', B0*v1*J1, quad_weights_2d, optimize='greedy')
    #             if isinstance(v1, (int, float)):
    #                 cewti2 = - np.einsum('ki, i -> k', B1*v0*J0, quad_weights_2d, optimize='greedy')
    #
    #             # noinspection PyUnboundLocalVariable
    #             GV[0, dofs] += np.concatenate((cewti1, cewti2))
    #
    #         else:
    #             raise NotImplementedError(f"can not handle non-orthogonal elements yet!")
    #
    #     GV = GlobalVector(GV.tocsr().T)
    #     return GV

    def ___PRIVATE_operator_wedge___(self, f2, quad_degree=None):
        """"""
        assert self.ndim == f2.ndim and self.k + f2.k == 3, " <___STORAGE_OPERATORS_WEDGE___> "
        assert self.mesh == f2.mesh, "Meshes do not match."

        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], f2.dqp[i]])) for i in
                           range(3)]
        quad_nodes, _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(
            quad_degree)
        _, bfSelf = self.do.evaluate_basis_at_meshgrid(
            *quad_nodes)
        _, bfOther = f2.do.evaluate_basis_at_meshgrid(
            *quad_nodes)

        W00 = self.___PRIVATE_wedge_H1___(quad_weights, bfOther[0], bfSelf[0])
        W11 = self.___PRIVATE_wedge_H1___(quad_weights, bfOther[1], bfSelf[1])
        W22 = self.___PRIVATE_wedge_H1___(quad_weights, bfOther[2], bfSelf[2])
        W = bmat([(W00, None, None),
                  (None, W11, None),
                  (None, None, W22)], format='csc')
        return W


    @staticmethod
    def ___PRIVATE_wedge_H1___(quad_weights, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights, bfO, bfS, optimize='optimal')
        return csc_matrix(M)


