# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly
from root.config.main import *


class _3dCSCG_Standard_Form_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()


    def RESET_cache(self):
        self._sf_.RESET_cache()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._sf_.___PRIVATE_do_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def evaluate_basis_at_quadrature(self, quad_degree, quad_type=None, compute_xietasigma=True):
        return self._sf_.space.do.evaluate_form_basis_at_quadrature(
            self._sf_.k, quad_degree, quad_type=quad_type, compute_xietasigma=compute_xietasigma)

    def compute_L2_energy_with(self, other=None, M=None):
        """Compute (self, other)_{L2} = int_{Omega}(self dot other)

        :param other: When it is None, other=self. Then we are doing self.compute_Ln_energy(n=2)
        :param M:
        :return:
        """
        if other is None: other = self._sf_
        assert self._sf_.mesh == other.mesh, "Meshes do not match."
        if M is None: M = self._sf_.operators.inner(other)

        if len(self._sf_.mesh.elements) == 0:
            LOCAL = 0

        elif self._sf_.mesh.elements.IS.homogeneous_according_to_types_wrt_metric:
            i = self._sf_.mesh.elements.indices[0]
            repM = M[i].toarray() # representative Mass matrix
            LOCAL = np.einsum('ij, ki, kj -> ',
                              repM,
                              self._sf_.cochain.array,
                              other.cochain.array,
                              optimize='greedy')
        else:
            LOCAL = list()
            for i in self._sf_.mesh.elements:
                LOCAL.append(self._sf_.cochain.local[i] @ M[i] @ other.cochain.local[i])
            LOCAL = np.sum(LOCAL)

        return COMM.allreduce(LOCAL, op=MPI.SUM)

    def compute_Ln_energy(self, n=2, quad_degree=None, vectorized=False):
        """So compute int_{Omega}( self ** n). When n = 2, it is equal to `self.compute_L2_energy_with()`.

        :param n:
        :param quad_degree:
        :param vectorized: {bool} Do we use vectorized reconstruction?
        :return:
        """
        sf = self._sf_
        assert n >= 1 and isinstance(n ,int), f'n={n} is wrong.'

        if quad_degree is None: quad_degree = [sf.dqp[i] + 1 for i in range(3)]
        quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

        if sf.k in (0, 3):
            OneOrThree = 1
        else:
            OneOrThree = 3

        #-------- vectorized --------------------------------------------------------
        if vectorized:  # we use the vectorized reconstruction
            reconstruction = sf.reconstruct(*quad_nodes, ravel=True, vectorized=True, value_only=True)

            if reconstruction[0] is None: # no local mesh elements

                total_energy = 0 # total energy in this core is 0

            else:
                reconstruction = np.sum([reconstruction[_]**n for _ in range(OneOrThree)], axis=0)
                xi, et, sg = np.meshgrid(*quad_nodes, indexing='ij')
                xi = xi.ravel('F')
                et = et.ravel('F')
                sg = sg.ravel('F')
                detJ = sf.mesh.elements.coordinate_transformation.vectorized.Jacobian(xi, et, sg)

                if sf.mesh.elements.IS.homogeneous_according_to_types_wrt_metric:

                    total_energy = np.einsum('kw, w, w -> ',
                                             reconstruction,
                                             detJ,
                                             quad_weights,
                                             optimize='optimal')

                else:

                    total_energy = np.einsum('kw, kw, w -> ',
                                             reconstruction,
                                             detJ,
                                             quad_weights,
                                             optimize='optimal')

            total_energy = COMM.allreduce(total_energy, op=MPI.SUM)

        #-------- non-vectorized --------------------------------------------------------
        else:

            Jacobian = sf.mesh.elements.coordinate_transformation.QUAD_1d.Jacobian(quad_degree, 'Gauss')

            _, SV = sf.reconstruct(*quad_nodes, ravel=True)
            local_energy = 0
            for i in sf.mesh.elements.indices:
                detJ = Jacobian[i]
                LEIntermediate = np.sum([SV[i][m]**n for m in range(OneOrThree)], axis=0)
                local_energy += np.sum( LEIntermediate * detJ * quad_weights )

            total_energy = COMM.allreduce(local_energy, op=MPI.SUM)

        return total_energy

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

    def compute_Ln_diff_from(self, other, n=2, quad_degree=None):
        """Compute: ||self - other||_{L^n} = n-root{ int_{Omega}(self - other)^n }.
        """
        sf = self._sf_
        of = other
        assert '3dCSCG_standard_form' in of.standard_properties.tags, "Other should be a _3dCSCG standard form."

        if sf.k in (0, 3):
            assert of.k in (0, 3)
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
            core_local = COMM.gather(core_local, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                globalError = np.sum(core_local) ** (1/n)
            else:
                globalError = None
            globalError = COMM.bcast(globalError, root=MASTER_RANK)
            return globalError
        else:
            raise NotImplementedError('Can only work on forms from the same mesh now.')

    def discretize(self, *args, **kwargs):
        return self._sf_.discretize(*args, **kwargs)

    def reconstruct(self, *args, **kwargs):
        return self._sf_.reconstruct(*args, **kwargs)

    def make_reconstruction_matrix_on_grid(self, xi, eta, sigma, element_range=None):
        """Make the reconstruction matrices for all mesh elements.

        These matrices are stored in a dict whose keys are the numbers of mesh elements and values
        are the local reconstruction matrices.

        Let `RM` be the reconstruction matrix (or the tuple of three matrices).
        If we want to do the local reconstruction, we do

            RM[i] @ f.cochain.local[i]

        and we will get the reconstructions of the form `f` on `meshgrid(xi, eta, sigma)` in mesh-element
        #i. And if `f` is a scalar form, we get a 1d array. And if `f` is a vector form, we get a
        tuple of three 1d arrays (its three components along x, y, z directions.)

        :param xi: 1d array
        :param eta: 1d array
        :param sigma: 1d array
        :param element_range:
            We are going to construct matrices for these mesh elements. It can be one of
                1) None: for all local elements
                2) 'mesh boundary': those local elements are attached to mesh boundary.

        """
        return self._sf_.___PRIVATE_make_reconstruction_matrix_on_grid___(xi, eta, sigma, element_range=element_range)

    @property
    def boundary_integrate(self):
        """"""
        return self._sf_._BI_