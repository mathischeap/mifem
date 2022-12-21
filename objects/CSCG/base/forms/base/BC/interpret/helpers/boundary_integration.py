# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/19/2022 4:22 PM
"""
from scipy.sparse import csc_matrix, csr_matrix
from components.freeze.main import FrozenOnly
import numpy as np
from components.quadrature import Quadrature
from components.assemblers import VectorAssembler


class CSCG_FORM_BC_Interpret_BoundaryIntegration(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._mesh_ = f.mesh
        self._elements_ = f.mesh.elements
        self.___empty___ = csc_matrix((f.num.basis, 1))
        self.quad_nodes, quad_weights = Quadrature(
            [i + 1 for i in self._f_.dqp], category='Gauss'
        ).quad

        x, y = np.meshgrid(self.quad_nodes[1], self.quad_nodes[2], indexing='ij')
        xNS = x.ravel('F')
        yNS = y.ravel('F')
        x, y = np.meshgrid(self.quad_nodes[0], self.quad_nodes[2], indexing='ij')
        xWE = x.ravel('F')
        yWE = y.ravel('F')
        x, y = np.meshgrid(self.quad_nodes[0], self.quad_nodes[1], indexing='ij')
        xBF = x.ravel('F')
        yBF = y.ravel('F')
        self._QN_ = {
            'N': [xNS, yNS],
            'S': [xNS, yNS],
            'W': [xWE, yWE],
            'E': [xWE, yWE],
            'B': [xBF, yBF],
            'F': [xBF, yBF],
        }

        self._qw_ = {
            'N': np.kron(quad_weights[2], quad_weights[1]),
            'S': np.kron(quad_weights[2], quad_weights[1]),
            'W': np.kron(quad_weights[2], quad_weights[0]),
            'E': np.kron(quad_weights[2], quad_weights[0]),
            'B': np.kron(quad_weights[1], quad_weights[0]),
            'F': np.kron(quad_weights[1], quad_weights[0]),
        }

        self.basis = self._f_.do.evaluate_basis_at_meshgrid(
            *self.quad_nodes, compute_xietasigma=False
        )[1]
        self._representing_boundaries_ = -1
        self._rm_ = self.___Pr_renew_RM___()
        if '3dCSCG_local_trace_form' in self._f_.standard_properties.tags:
            self._assembler_ = VectorAssembler(self._f_.numbering.local_gathering)
        else:
            pass
        self._freeze_self_()

    def ___Pr_renew_RM___(self):
        """"""
        if self._representing_boundaries_ == id(self._f_.BC.boundaries):
            return self._rm_   # no need to renew the reconstruction matrices

        else:   # renew rm

            rm = self._f_.reconstruct.___PrLT_make_reconstruction_matrix_on_grid___(
                *self.quad_nodes, element_range=self._f_.BC._involved_element_parts_
            )

            self._representing_boundaries_ = id(self._f_.BC.boundaries)
            self._rm_ = rm
            return rm

    def __call__(self, i):
        """Return the boundary local cochains in real time (following the current `BC.CF` and `BC.boundaries`).
        for mesh-element i. For those dofs not locating on the mesh boundary, we set 0 for them.

        Parameters
        ----------
        i

        Returns
        -------

        """
        if self._f_.BC.boundaries == list() or \
           self._f_.BC.boundaries is None:

            return self.___empty___

        else:

            rm = self.___Pr_renew_RM___()

            if i not in rm:
                return self.___empty___

            else:

                if '3dCSCG_local_trace_form' in self._f_.standard_properties.tags:

                    cf = self._f_.BC.CF
                    valid_sides = rm[i].keys()
                    Tmap = self._mesh_.trace.elements.map[i]

                    local_vector = dict()

                    if cf.ftype == 'boundary-wise':

                        for vs in 'NSWEBF':
                            if vs in valid_sides:
                                e = Tmap['NSWEBF'.index(vs)]
                                te = self._mesh_.trace.elements[e]
                                x, y, z = te.coordinate_transformation.mapping(*self._QN_[vs])
                                boundary_name = te.on_mesh_boundary
                                func = self._f_.BC.CF.do.evaluate_func_at_time()[boundary_name]
                                func = [_(x, y, z) for _ in func]
                                basis = self.basis[vs]
                                constant_Jacobian = te.coordinate_transformation.constant.Jacobian

                                if len(func) == 1:  # we are looking at a scalar.
                                    func = func[0]
                                    basis = basis[0]

                                    if constant_Jacobian is None:
                                        raise NotImplementedError()

                                    else:
                                        element_side_vector = np.einsum(
                                            'ij, j -> i',
                                            basis,
                                            func * self._qw_[vs] * constant_Jacobian,
                                            optimize='greedy',
                                        )

                                        local_vector[vs] = element_side_vector

                                else:
                                    raise NotImplementedError()

                            else:
                                local_vector[vs] = np.zeros(self._f_.num.basis_onside[vs])

                        spa_vec = self._assembler_(local_vector, 'add')

                        spa_vec = csr_matrix(spa_vec).T

                        return spa_vec

                    else:
                        raise NotImplementedError()

                else:
                    raise NotImplementedError()

        #     element = self._elements_[i]
        #     if element.whether.internal:
        #         # cause this is a realtime function, we always check self._f_.BC.
        #         return self.___empty___
        #
        #     # else:
        #     #     if i in self.___Pr_BcCo___:
        #     #         return csr_matrix(self.___Pr_BcCo___[i]).T
        #     #     else:
        #     #         return self.___empty___
