import numpy as np

from screws.freeze.inheriting.frozen_only import FrozenOnly
from root.config.main import cOmm, MPI


class _2dCSCG_Standard_Form_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def reset_cache(self):
        self._sf_.___PRIVATE_reset_cache___()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._sf_.___PRIVATE_do_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def discretize(self, *args, **kwargs):
        return self._sf_.discretize(*args, **kwargs)

    def reconstruct(self, *args, **kwargs):
        return self._sf_.reconstruct(*args, **kwargs)

    def make_reconstruction_matrix_on_grid(self, xi, eta):
        """Make the reconstruction matrices for all mesh elements. These matrices are stored in
        a dict whose keys are the numbers of mesh elements and values are the local reconstruction
        matrices.

        Let `RM` be the reconstruction matrix. If we want to do the local reconstruction, we do

            RM[i] @ f.cochain.local[i]

        and we will get the reconstructions of the form `f` on `meshgrid(xi, eta)` in mesh-element
        #i. And if `f` is a scalar form, we get a 1d array. And if `f` is a vector form, we get a
        tuple of two 1d arrays (its two components along x, y directions.)

        """
        return self._sf_.___PRIVATE_make_reconstruction_matrix_on_grid___(xi, eta)

    def compute_Ln_energy(self, n=2, quad_degree=None):
        """Compute int_{Omega}(self^n).

        :param n: default n=2, we compute the L2 inner product between self and self: (self, self)_{L2}.
        :param quad_degree:
        :return:
        """
        f = self._sf_
        if quad_degree is None:
            if n <= 1:
                quad_degree = f.dqp
            elif n == 'infinity':
                quad_degree = [int(f.dqp[_] * 2) for _ in range(2)]
            else:
                quad_degree = [int(f.dqp[_] * (1 + (n-1)/2)) for _ in range(2)]
        else:
            pass
        assert len(quad_degree) == 2 and all([isinstance(_, int) for _ in quad_degree]), \
            f"quad_degree = {quad_degree} is illegal."

        quad_nodes, _, quad_weights_1d = \
            f.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        xi, et = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel()
        et = et.ravel()
        reconstruction = f.reconstruct(*quad_nodes, ravel=True)[1]
        detJ = f.mesh.elements.coordinate_transformation.Jacobian(xi, et)

        total_energy = 0
        if f.k in (0, 2): # scalar form
            for i in f.mesh.elements:

                if n == 'infinity':
                    Ei = np.max(np.abs(reconstruction[i][0]))

                else:
                    v = reconstruction[i][0]**n
                    Ei = np.sum(v * quad_weights_1d * detJ[i])

                total_energy += Ei
        else:
            raise NotImplementedError('Not implemented for 1-form.')

        total_energy = cOmm.allreduce(total_energy, op=MPI.SUM)

        return total_energy

    def compute_L2_energy_with(self, other=None, M=None):
        """Compute (self, other)_{L2} = int_{Omega}(self dot other)

        :param other: When it is None, other=self. Then we are doing self.compute_Ln_energy(n=2)
        :param M:
        :return:
        """
        if other is None: other = self._sf_
        assert self._sf_.mesh == other.mesh, "Meshes do not match."
        if M is None: M = self._sf_.operators.inner(other)
        LOCAL = list()
        for i in self._sf_.mesh.elements:
            LOCAL.append(self._sf_.cochain.local[i] @ M[i] @ other.cochain.local[i])
        LOCAL = np.sum(LOCAL)
        return cOmm.allreduce(LOCAL, op=MPI.SUM)

    def compute_Ln_norm(self, n=2, quad_degree=None):
        """Compute  ||self ||_{L^n} = ( int_{Omega}(self^n) )^(1/n) , which is  the n-root of Ln-energy.

        :param n: {default n=2} ||self ||_{L^2}
        :param quad_degree:
        :return:
        """
        total_energy = self.compute_Ln_energy(n=n, quad_degree=quad_degree)
        return total_energy**(1/n)

    def compute_Ln_norm_of_coboundary(self, n=2, quad_degree=None):
        """Compute || d(self) ||_{L^n} .

        We can, for example, use this method to conveniently evaluate how good the conservation of
        mass  d · u = 0 is satisfied with following measures:

            || d · u ||_{L^2} or || d · u ||_{L^infinity}

        """
        d_self = self._sf_.coboundary()

        Ln_norm_of_d_self = d_self.do.compute_Ln_norm(n=n, quad_degree=quad_degree)

        return Ln_norm_of_d_self