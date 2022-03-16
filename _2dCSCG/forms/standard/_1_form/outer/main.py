# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from _2dCSCG.forms.standard._1_form.outer.special import _1Form_Outer_Special
import numpy as np
from scipy import sparse as spspa


from _2dCSCG.forms.standard._1_form.base.main import _1Form_BASE
from _2dCSCG.forms.standard._1_form.outer.discretize.main import _2dCSCG_S1Fo_Discretize



class _2dCSCG_1Form_Outer(_1Form_BASE):
    """
    Standard outer 1-form.

    :param mesh:
    :param space:
    :param is_hybrid:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, is_hybrid=True,
        numbering_parameters='Naive',  name='outer-oriented-1-form'):
        super().__init__(mesh, space, is_hybrid, 'outer', numbering_parameters, name)
        super().__init_1form_base__()
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_1form_Outer')
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_1form')
        self._special_ = _1Form_Outer_Special(self)
        self._discretize_ = _2dCSCG_S1Fo_Discretize(self)
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    @property
    def special(self):
        return self._special_

    @property
    def discretize(self):
        return self._discretize_

    def reconstruct(self, xi, eta, ravel=False, i=None):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param i:
        :return:
        """
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, eta)
        xyz = dict()
        value = dict()
        shape = [len(xi), len(eta)]
        INDICES = self.mesh.elements.indices if i is None else [i,]
        iJ = self.mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
        for i in INDICES:
            element = self.mesh.elements[i]
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            u = np.einsum('ij, i -> j', basis[0], self.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
            v = np.einsum('ij, i -> j', basis[1], self.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')
            value[i] = [None, None]
            typeWr2Metric = element.type_wrt_metric.mark
            iJi = iJ[i]
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                value[i][0] = u*iJi[1][1]
                value[i][1] = v*iJi[0][0]
            else:
                value[i][0] = + u*iJi[1][1] - v*iJi[0][1]
                value[i][1] = - u*iJi[1][0] + v*iJi[0][0]
            if ravel:
                pass
            else:
                # noinspection PyUnresolvedReferences
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                # noinspection PyUnresolvedReferences
                value[i] = [value[i][j].reshape(shape, order='F') for j in range(2)]
        return xyz, value

    def ___PRIVATE_make_reconstruction_matrix_on_grid___(self, xi, eta):
        """
        Make a dict (keys are #mesh-elements) of matrices whose columns refer to
        nodes of meshgrid(xi, eta, indexing='ij') and rows refer to
        local dofs.

        If we apply these matrices to the local dofs, we will get the
        reconstructions on the nodes in the mesh-elements.

        :param xi: 1d array in [-1, 1].
        :param eta: 1d array in [-1, 1].
        :return:
        """
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, eta)

        INDICES = self.mesh.elements.indices
        iJ = self.mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)

        b0, b1 = basis[0].T, basis[1].T
        OO01 = 0 * b1
        OO10 = 0 * b0

        type_cache = dict()
        RM = dict()
        for i in INDICES:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in type_cache:
                    RM[i] = type_cache[typeWr2Metric]
                else:
                    iJi = iJ[i]
                    rm00 = + np.einsum('ji, j -> ji', b0, iJi[1][1], optimize='greedy')
                    rm11 = + np.einsum('ji, j -> ji', b1, iJi[0][0], optimize='greedy')
                    if typeWr2Metric[:4] == 'Orth':
                        RM_i_ = ( np.hstack((rm00, OO01)),
                                  np.hstack((OO10, rm11)) )
                    else:
                        rm01 = - np.einsum('ji, j -> ji', b1, iJi[0][1], optimize='greedy')
                        rm10 = - np.einsum('ji, j -> ji', b0, iJi[1][0], optimize='greedy')
                        RM_i_ = ( np.hstack((rm00, rm01)),
                                  np.hstack((rm10, rm11)) )

                    type_cache[typeWr2Metric] = RM_i_
                    RM[i] = RM_i_

            else:
                iJi = iJ[i]
                rm00 = + np.einsum('ji, j -> ji', b0, iJi[1][1], optimize='greedy')
                rm01 = - np.einsum('ji, j -> ji', b1, iJi[0][1], optimize='greedy')
                rm10 = - np.einsum('ji, j -> ji', b0, iJi[1][0], optimize='greedy')
                rm11 = + np.einsum('ji, j -> ji', b1, iJi[0][0], optimize='greedy')
                RM[i] = ( np.hstack((rm00, rm01)),
                          np.hstack((rm10, rm11)) )

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
        M00 = self.___PRIVATE_inner_Helper1___(quad_weights, sqrtg*g[1][1], bfOther[0], bfSelf[0])
        M11 = self.___PRIVATE_inner_Helper1___(quad_weights, sqrtg*g[0][0], bfOther[1], bfSelf[1])
        if isinstance(mark, str) and mark[:4] == 'Orth':
            M01 = None
            M10 = None
        else:
            M01 = self.___PRIVATE_inner_Helper1___(quad_weights, -sqrtg*g[1][0], bfOther[0], bfSelf[1])
            if other is self:
                M10 = M01.T
            else:
                M10 = self.___PRIVATE_inner_Helper1___(quad_weights, -sqrtg*g[0][1], bfOther[1], bfSelf[0])
        Mi = spspa.bmat([(M00, M01),
                         (M10, M11)], format='csc')
        return Mi

    @staticmethod
    def ___PRIVATE_inner_Helper1___(quad_weights, sqrtg_g, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights*sqrtg_g, bfO, bfS, optimize='greedy')
        return spspa.csc_matrix(M)


    def ___PRIVATE_operator_wedge___(self, other, quad_degree=None):
        """ """
        assert other.k == 1, "Need a _2dCSCG_1Form"
        assert self.mesh == other.mesh, "Meshes do not match."
        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], other.dqp[i]])) for i in range(2)]
        quad_nodes, _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        xietasigma, bS = self.do.evaluate_basis_at_meshgrid(*quad_nodes)
        _, bO = other.do.evaluate_basis_at_meshgrid(*quad_nodes)
        W00 = np.einsum('im, jm -> ij', bO[0], bS[0]*quad_weights[np.newaxis, :], optimize='optimal')
        W11 = np.einsum('im, jm -> ij', bO[1], bS[1]*quad_weights[np.newaxis, :], optimize='optimal')
        i, j = other.num.basis_components
        m, n = self.num.basis_components
        #      m   n
        # i  |W00 W01 |
        # j  |W10 W11 |
        W = np.vstack((np.hstack((W00, np.zeros((i,n)))),
                       np.hstack((np.zeros((j,m)), W11))))
        return spspa.csc_matrix(W)




if __name__ == '__main__':
    # mpiexec python _2dCSCG\forms\standard\_1_form\outer\main.py
    from _2dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    # mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1],[0,1]))([10,10])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',2)])
    FC = FormCaller(mesh, space)

    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    f1 = FC('1-f-o', is_hybrid=True)

    # M0 = f1.matrices.mass[0]
    # print(M0.toarray())
    #
    # E21 = f1.matrices.incidence[0]
    # print(E21.toarray())

    f1.TW.func.do.set_func_body_as(ES, 'velocity')
    f1.TW.current_time = 0
    f1.TW.do.push_all_to_instant()
    f1.discretize()
    print(f1.error.L())
    #
    # from root.mifem import save
    #
    # save(f1, 'test_2d_f1_o')