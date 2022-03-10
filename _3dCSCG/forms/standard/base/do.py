
from screws.freeze.main import FrozenOnly
from root.config.main import *

class _3dCSCG_Standard_Form_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def reset_cache(self):
        self._sf_.___PRIVATE_reset_cache___()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._sf_.___PRIVATE_do_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def evaluate_basis_at_quadrature(self, quad_degree, quad_type=None, compute_xietasigma=True):
        return self._sf_.space.do.evaluate_form_basis_at_quadrature(
            self._sf_.k, quad_degree, quad_type=quad_type, compute_xietasigma=compute_xietasigma)

    def resemble(self, *args, **kwargs):
        return self._sf_.___PRIVATE_do_resemble___(*args, **kwargs)

    def interpolate(self, data_set):
        """Like resemble, but `interpolate` is from a data set, not another form.

        :param data_set:
        :return:
        """
        raise NotImplementedError()

    def cross_product(self, *args, **kwargs):
        return self._sf_.special.cross_product(*args, **kwargs)

    def compute_L2_energy_with(self, other=None, M=None):
        """Compute (self, other)_{L2} = int_{Omega}(self dot other)

        :param other: When it is None, other=self. Then we are doing self.compute_Ln_energy(n=2)
        :param M:
        :return:
        """
        return self._sf_.___PRIVATE_do_compute_Ln_energy_with___(other=other, M=M)

    def compute_Ln_diff_from(self, other, n=2, quad_degree=None):
        """ compute: ||self - other||_{L^n} = n-root{ int_{Omega}(self - other)^n }."""
        sf = self._sf_
        of = other
        assert '3dCSCG_standard_form' in of.standard_properties.tags, "Other should be a _3dCSCG standard form."
        if sf.k in (0, 3):
            assert sf.k in (0, 3)
            OneOrThree = 1
        else:
            assert of.k in (1, 2)
            OneOrThree = 3
        if sf.mesh == of.mesh:
            if quad_degree is None: quad_degree = [np.max([sf.dqp[i], of.dqp[i]]) + 1 for i in range(3)]
            quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
            Jacobian = sf.mesh.elements.coordinate_transformation.QUAD_1d.Jacobian(quad_degree, 'Gauss')
            _, SV = sf.reconstruct(*quad_nodes, ravel=True)
            _, OV = of.reconstruct(*quad_nodes, ravel=True)
            localError = list()
            for i in sf.mesh.elements.indices:
                detJ = Jacobian[i]
                LEIntermediate = np.sum([(SV[i][m] - OV[i][m])**n for m in range(OneOrThree)], axis=0)
                localError.append(np.sum(LEIntermediate * detJ * quad_weights))
            core_local = np.sum(localError)
            core_local = cOmm.gather(core_local, root=mAster_rank)
            if rAnk == mAster_rank:
                globalError = np.sum(core_local) ** (1/n)
            else:
                globalError = None
            globalError = cOmm.bcast(globalError, root=mAster_rank)
            return globalError
        else:
            raise NotImplementedError('Can only work on forms from the same mesh now.')

    def discretize(self, *args, **kwargs):
        return self._sf_.discretize(*args, **kwargs)

    def reconstruct(self, *args, **kwargs):
        return self._sf_.reconstruct(*args, **kwargs)

    def make_reconstruction_matrix_on_grid(self, xi, eta, sigma):
        """Make the reconstruction matrices for all mesh elements. These matrices are stored in
        a dict whose keys are the numbers of mesh elements and values are the local reconstruction
        matrices.

        Let `RM` be the reconstruction matrix (or the tuple of three matrices).
        If we want to do the local reconstruction, we do

            RM[i] @ f.cochain.local[i]

        and we will get the reconstructions of the form `f` on `meshgrid(xi, eta, sigma)` in mesh-element
        #i. And if `f` is a scalar form, we get a 1d array. And if `f` is a vector form, we get a
        tuple of three 1d arrays (its three components along x, y, z directions.)

        """
        return self._sf_.___PRIVATE_make_reconstruction_matrix_on_grid___(xi, eta, sigma)