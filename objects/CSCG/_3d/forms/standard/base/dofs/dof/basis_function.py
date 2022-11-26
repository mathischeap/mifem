# -*- coding: utf-8 -*-
"""The class for the basis function of a dof (not dofs)."""


from components.freeze.main import FrozenOnly
import numpy as np


class _3dCSCG_SF_DOF_BF(FrozenOnly):
    """"""
    def __init__(self, dof):
        self._dof_ = dof
        self._sf_ = dof._sf_
        self._mesh_ = self._sf_.mesh
        self._space_ = self._sf_.space
        self._freeze_self_()


    def reconstruct(self, xi, et, sg, ravel=False):
        """We reconstruct this single basis function in the corresponding local mesh element(s).

        :param xi: must be increasing 1d array in [-1,1].
        :param et: must be increasing 1d array in [-1,1].
        :param sg: must be increasing 1d array in [-1,1].
        :param ravel: If ravelled, we return flat results, otherwise, we return 3d results.
        :return: A tuple of reconstructions corresponding to the positions.
        """
        positions = self._dof_.positions # get all local positions.

        if positions == list():
            return tuple(), tuple()

        k = self._sf_.k
        mesh = self._mesh_
        space = self._space_
        XYZ = tuple()
        IN_SITE_BF = tuple()

        shape = [len(xi), len(et), len(sg)]
        xietasigma, RBF= space.DO_evaluate_form_basis_at_meshgrid(k, xi, et ,sg, compute_xietasigma=True)
        RBF = np.vstack(RBF)
        for POS in positions:
            E, I = POS
            rbf = RBF[I] # find the bf for this position
            ME = mesh.elements[E] # the mesh-element of this position.

            #----------- do the reconstruction according to the k ------------------------------------------
            if k == 0:
                xyz, in_site_bf = self.___PRIVATE_reconstruct_0_bf___(xietasigma, rbf, ME)
            else:
                raise NotImplementedError()
            #---------- post process the data ---------------------------------------------------------------
            if ravel:
                pass
            else:
                xyz = [xyz[_].reshape(shape, order='F') for _ in range(3)]
                in_site_bf = [in_site_bf[_].reshape(shape, order='F') for _ in range(len(in_site_bf))]
            #============================================================================================

            XYZ += (xyz,)
            IN_SITE_BF += (in_site_bf,)

        return XYZ, IN_SITE_BF

    @staticmethod
    def ___PRIVATE_reconstruct_0_bf___(xi_et_sg, rbf, ME):
        """We make the in-site bf from the rbf (reference bf) for the 0-standard-form.

        :param xi_et_sg: rbf is evaluated at meshgrid(*xi_et_sg).
        :param rbf: reference bf
        :param ME: this bf is in this Mesh-Element #.
        :return: The in-site bf: xyz and bfv (basis function value).
        """
        xyz = ME.coordinate_transformation.mapping(*xi_et_sg)
        return xyz, (rbf,) # for 0-standard-form, in-site-bf is equal to reference bf.
    def ___PRIVATE_reconstruct_1_bf___(self, xi_et_sg, rbf, ME):
        """We make the in-site bf from the rbf (reference bf) for the 1-standard-form.

        :param xi_et_sg: rbf is evaluated at meshgrid(*xi_et_sg).
        :param rbf: reference bf
        :param ME: this bf is in this Mesh-Element #.
        :return: The in-site bf: xyz and bfv (basis function value).
        """
        raise NotImplementedError()
    def ___PRIVATE_reconstruct_2_bf___(self, xi_et_sg, rbf, ME):
        """We make the in-site bf from the rbf (reference bf) for the 2-standard-form.

        :param xi_et_sg: rbf is evaluated at meshgrid(*xi_et_sg).
        :param rbf: reference bf
        :param ME: this bf is in this Mesh-Element #.
        :return: The in-site bf: xyz and bfv (basis function value).
        """
        raise NotImplementedError()
    def ___PRIVATE_reconstruct_3_bf___(self, xi_et_sg, rbf, ME):
        """We make the in-site bf from the rbf (reference bf) for the 3-standard-form.

        :param xi_et_sg: rbf is evaluated at meshgrid(*xi_et_sg).
        :param rbf: reference bf
        :param ME: this bf is in this Mesh-Element #.
        :return: The in-site bf: xyz and bfv (basis function value).
        """
        raise NotImplementedError()




    def visualize(self, *args, **kwargs):
        """We use matplotlib to visualize the basis function of this dof in-site.

        "in-site" means we will visualize it in the corresponding mesh element instead of the
        reference domain.

        Cause a dof instance is already local (may in multiple cores), we will call this visualizing
        method locally.
        """
        return getattr(self, f"___PRIVATE_visualize_{self._sf_.k}form___")(*args, **kwargs)

    def ___PRIVATE_visualize_0form___(self, density=1000):
        """"""
        # density = int(density**(1/3)) + 1
        # xi = eta = sigma = np.linspace(-1, 1, density)
        raise NotImplementedError()
    def ___PRIVATE_visualize_1form___(self, density=1000):
        """"""
        raise NotImplementedError()
    def ___PRIVATE_visualize_2form___(self, density=1000):
        """"""
        raise NotImplementedError()
    def ___PRIVATE_visualize_3form___(self, density=1000):
        """"""
        raise NotImplementedError()
