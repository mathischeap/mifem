# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:34 PM
"""
import sys

import numpy as np

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from root.config.main import MPI, COMM


class miUs_Triangular_SF_Do(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        """"""
        xi_eta, basis = getattr(self._sf_.space.evaluation,
                                self._sf_.__class__.__name__)(*args, **kwargs)

        if self._sf_.k in (0, 2):
            return xi_eta, basis

        else:
            BASIS_dict = dict()
            tM = self._sf_.IDT.transition_matrices
            tT = self._sf_.IDT.transition_types
            b0, b1 = basis
            num0 = b0.shape[0]
            BASIS = np.vstack(basis)

            basis_pool = dict()

            for e in self._sf_.mesh.elements.indices:
                if e not in tM:
                    BASIS_dict[e] = basis
                else:
                    key = tT[e]
                    if key in basis_pool:
                        b_ = basis_pool[key]
                    else:
                        BASIS_d = tM[e] @ BASIS
                        b_ = [BASIS_d[:num0, :], BASIS_d[num0:]]
                        basis_pool[key] = b_

                    BASIS_dict[e] = b_

            return xi_eta, BASIS_dict


    def make_reconstruction_matrix_on_grid(self, xi, eta, element_range=None):
        """Make the reconstruction matrices for all mesh elements. These matrices are stored in
        a dict whose keys are the numbers of mesh elements and values are the local reconstruction
        matrices.

        Let `RM` be the reconstruction matrix. If we want to do the local reconstruction, we do

            RM[i] @ f.cochain.local[i]

        and we will get the reconstructions of the form `f` on `meshgrid(xi, eta)` in mesh-element
        #i. And if `f` is a scalar form, we get a 1d array. And if `f` is a vector form, we get a
        tuple of two 1d arrays (its two components along x, y directions.)

        """

        if element_range is None:
            INDICES = self._sf_.mesh.elements.indices
        elif element_range.__class__.__name__ in ('int', 'float', 'int32', 'int64'):
            INDICES = [int(element_range),]
        else:
            raise Exception()

        xietasigma, basis = self.evaluate_basis_at_meshgrid(xi, eta)
        elements = self._sf_.mesh.elements
        RM = dict()

        if self._sf_.k == 0:
            for i in INDICES:
                RM[i] = basis[0].T

        elif self._sf_.k == 1:
            if self._sf_.orientation == 'inner':
                for i in INDICES:
                    element = elements[i]
                    b0, b1 = basis[i][0].T, basis[i][1].T
                    iJi = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
                    rm00 = np.einsum('ji, j -> ji', b0, iJi[0][0], optimize='optimal')
                    rm01 = np.einsum('ji, j -> ji', b1, iJi[1][0], optimize='optimal')
                    rm10 = np.einsum('ji, j -> ji', b0, iJi[0][1], optimize='optimal')
                    rm11 = np.einsum('ji, j -> ji', b1, iJi[1][1], optimize='optimal')
                    RM[i] = (np.hstack((rm00, rm01)),
                             np.hstack((rm10, rm11)))

            else:
                for i in INDICES:
                    element = elements[i]
                    b0, b1 = basis[i][0].T, basis[i][1].T
                    iJi = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
                    rm00 = + np.einsum('ji, j -> ji', b0, iJi[1][1], optimize='greedy')
                    rm01 = - np.einsum('ji, j -> ji', b1, iJi[0][1], optimize='greedy')
                    rm10 = - np.einsum('ji, j -> ji', b0, iJi[1][0], optimize='greedy')
                    rm11 = + np.einsum('ji, j -> ji', b1, iJi[0][0], optimize='greedy')
                    RM[i] = (np.hstack((rm00, rm01)),
                             np.hstack((rm10, rm11)))

        elif self._sf_.k == 2:
            for i in INDICES:
                element = elements[i]
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                RM[i] = np.einsum('ij, j -> ji', basis[0], det_iJ, optimize='greedy')

        else:
            raise Exception()

        return RM



    def compute_Ln_energy(self, n=2, quad_degree=None):
        """Compute int_{Omega}(self^n).

        :param n: default n=2, we compute the L2 inner product between self and self: (self, self)_{L2}.
        :param quad_degree:
        :return:
        """
        f = self._sf_
        assert n >= 1 and isinstance(n ,int), f'n={n} is wrong.'
        #-------- parse quad_degree ---------------------------------------------------
        if quad_degree is None:
            quad_degree = [f.space.p + 1 for _ in range(2)]
        else:
            pass

        quad_nodes, _, quad_weights_1d = \
            f.space.evaluation.quadrature(quad_degree)

        xi, et = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        et = et.ravel('F')
        reconstruction = f.reconstruct(*quad_nodes, ravel=True)[1]

        total_energy = 0
        if f.k in (0, 2): # scalar form
            for i in f.mesh.elements:
                detJ = f.mesh.elements[i].coordinate_transformation.Jacobian(xi, et)
                v = reconstruction[i][0]**n
                Ei = np.sum(v * quad_weights_1d * detJ)
                total_energy += Ei

        else:  # vector form
            for i in f.mesh.elements:
                detJ = f.mesh.elements[i].coordinate_transformation.Jacobian(xi, et)
                v = reconstruction[i][0]**n + reconstruction[i][1]**n
                Ei = np.sum(v * quad_weights_1d * detJ)
                total_energy += Ei

        total_energy = COMM.allreduce(total_energy, op=MPI.SUM)

        return total_energy

    def compute_Ln_norm(self, n=2, quad_degree=None):
        """Compute  ||self ||_{L^n} = ( int_{Omega}(self^n) )^(1/n) , which is the n-root of Ln-energy.

        :param n: {default n=2} ||self ||_{L^2}
        :param quad_degree:
        :return:
        """
        total_energy = self.compute_Ln_energy(n=n, quad_degree=quad_degree)
        return total_energy**(1/n)


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

        else:
            LOCAL = list()
            for i in self._sf_.mesh.elements:
                LOCAL.append(self._sf_.cochain.local[i] @ M[i] @ other.cochain.local[i])
            LOCAL = np.sum(LOCAL)

        return COMM.allreduce(LOCAL, op=MPI.SUM)

    def compute_Ln_norm_of_coboundary(self, n=2, quad_degree=None):
        """Compute || d(self) ||_{L^n} .

        We can, for example, use this method to conveniently evaluate how good the conservation of
        mass  d · u = 0 is satisfied with following measures:

            || d · u ||_{L^2} or || d · u ||_{L^infinity}

        """
        d_self = self._sf_.coboundary()

        Ln_norm_of_d_self = d_self.do.compute_Ln_norm(n=n, quad_degree=quad_degree)

        return Ln_norm_of_d_self


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
