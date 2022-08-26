# -*- coding: utf-8 -*-

from scipy.sparse import csc_matrix
import numpy as np


# noinspection PyUnresolvedReferences
class _3dCSCG_S0F_Private:
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
        _, basis = self.do.evaluate_basis_at_meshgrid(xi, et, sg)
        RM = dict()

        #------ take care of element range -----------------------------------------------
        if element_range is None:
            INDICES = self.mesh.elements
        elif element_range == 'mesh boundary':
            INDICES = self.mesh.boundaries.involved_elements
        else:
            raise Exception(f"element_range = {element_range} is wrong!")

        #-------- generating the reconstruction matrices for included mesh elements -------------
        rmi = basis[0].T
        for i in INDICES:
            RM[i] = rmi
        return RM

    def ___PRIVATE_operator_inner___(self, _, i, xietasigma, quad_weights, bfSelf, bfOther):
        """Note that here we only return a local matrix."""
        element = self.mesh.elements[i]
        detJ = element.coordinate_transformation.Jacobian(*xietasigma)
        Mi = np.einsum('im, jm, m -> ij', bfOther[0], bfSelf[0], detJ*quad_weights, optimize='greedy')
        Mi = csc_matrix(Mi)
        return Mi
