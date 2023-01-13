# -*- coding: utf-8 -*-
from root.config.main import *
from scipy import sparse as spspa
from objects.CSCG._2d.forms.standard.base.main import _2dCSCG_Standard_Form
from objects.CSCG._2d.forms.standard._0_form.base.discretize.main import _2dCSCG_S0F_Discretize
from objects.CSCG._2d.forms.standard._0_form.base.visualize.main import _2dCSCG_S0F_VIS
from objects.CSCG._2d.forms.standard._0_form.base.reconstruct import _2dCSCG_S0F_Reconstruct


class _0Form_BASE(_2dCSCG_Standard_Form):
    """"""
    def __init_0form_base__(self):
        self._discretize_ = _2dCSCG_S0F_Discretize(self)
        self._visualize_ = _2dCSCG_S0F_VIS(self)
        self._reconstruct_ = None

    @property
    def visualize(self):
        return self._visualize_

    def ___Pr_check_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2

        if func_body.__class__.__name__ == '_2dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"2dCSCG 0form FUNC does not accept func _2dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"2dCSCG 0form FUNC does not accept func {func_body.__class__}")

    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2
        raise Exception(f"2dCSCG 0form BC does not accept func {func_body.__class__}")

    @property
    def discretize(self):
        return self._discretize_

    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            # noinspection PyAttributeOutsideInit
            self._reconstruct_ = _2dCSCG_S0F_Reconstruct(self)
        return self._reconstruct_

    def ___PRIVATE_make_reconstruction_matrix_on_grid___(self, xi, eta, element_range=None):
        """
        Make a dict (keys are #mesh-elements) of matrices whose columns refer to
        nodes of meshgrid(xi, eta, indexing='ij') and rows refer to
        local dofs.

        If we apply these matrices to the local dofs, we will get the
        reconstructions on the nodes in the mesh-elements.

        :param xi: 1d array in [-1, 1]
        :param eta: 1d array in [-1, 1]
        :return:
        """
        _, basis = self.do.evaluate_basis_at_meshgrid(xi, eta)
        RM = dict()
        if element_range is None:
            INDICES = self.mesh.elements.indices
        elif element_range == 'mesh boundary':
            INDICES = self.mesh.boundaries.involved_elements
        else:
            raise Exception(f"element_range = {element_range} is wrong!")

        rmi = basis[0].T
        for i in INDICES:
            RM[i] = rmi
        return RM

    def ___PRIVATE_operator_inner___(self, _, i, xietasigma, quad_weights, bfSelf, bfOther):
        """Note that here we only return a local matrix."""
        element = self.mesh.elements[i]
        detJ = element.coordinate_transformation.Jacobian(*xietasigma)
        Mi = np.einsum('im, jm, m -> ij', bfOther[0], bfSelf[0], detJ*quad_weights, optimize='greedy')
        Mi = spspa.csc_matrix(Mi)
        return Mi

    def ___PRIVATE_operator_wedge___(self, other, quad_degree=None):
        """ """
        assert self.ndim == other.ndim and self.k + other.k == other.ndim, " <_0Form> "
        try:
            assert self.mesh == other.mesh
        except AssertionError:
            raise Exception(' <_0Form_int_wedge> : meshes do not fit.')
        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], other.dqp[i]])) + 1 for i in range(2)]
        quad_nodes, _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        xietasigma, basisS = self.do.evaluate_basis_at_meshgrid(*quad_nodes)
        _, basisO = other.do.evaluate_basis_at_meshgrid(*quad_nodes)
        W = np.einsum('im, jm, m -> ij', basisO[0], basisS[0], quad_weights, optimize='greedy')
        return spspa.csc_matrix(W)
