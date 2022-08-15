# -*- coding: utf-8 -*-


from root.config.main import *
from scipy.sparse import csc_matrix


# noinspection PyUnresolvedReferences
class _3dCSCG_S3F_Private:

    def ___PRIVATE_make_reconstruction_matrix_on_grid___(self, xi, et, sg):
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
        :return:
        """
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, et, sg)
        biJC = dict()
        RM = dict()
        for i in self.mesh.elements:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in biJC:
                    RM[i] = biJC[typeWr2Metric]
                else:
                    det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                    basis_det_iJ = (basis[0] * det_iJ).T
                    biJC[typeWr2Metric] = basis_det_iJ
                    RM[i] = basis_det_iJ
            else:
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                basis_det_iJ = basis[0] * det_iJ
                RM[i] = basis_det_iJ.T
        return RM

    def ___PRIVATE_operator_inner___(self, _, i, xietasigma, quad_weights, bfSelf, bfOther):
        """Note that here we only return a local matrix."""
        element = self.mesh.elements[i]
        detJ = element.coordinate_transformation.Jacobian(*xietasigma)
        Mi = np.einsum('im, jm, m -> ij',
            bfOther[0], bfSelf[0], np.reciprocal(detJ)*quad_weights,
            optimize='greedy'
        )
        Mi = csc_matrix(Mi)
        return Mi

    def ___PRIVATE_operator_wedge___(self, other, quad_degree=None):
        """In fact, it is integral over wedge product."""
        assert other.__class__.__name__ == '_0Form', "Need a _3dCSCG_0Form"
        assert self.mesh == other.mesh, "Meshes do not match."
        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], other.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        _, basisS = self.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        _, basisO = other.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        W = np.einsum('im, jm, m -> ij', basisO[0], basisS[0], quad_weights, optimize='greedy')
        return csc_matrix(W)
