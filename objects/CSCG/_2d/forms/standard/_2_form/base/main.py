# -*- coding: utf-8 -*-
from root.config.main import *
from scipy import sparse as spspa
from objects.CSCG._2d.forms.standard.base.main import _2dCSCG_Standard_Form
from objects.CSCG._2d.forms.standard._2_form.base.discretize.main import _2dCSCG_S2F_Discretize
from objects.CSCG._2d.forms.standard._2_form.base.visualize.main import _2dCSCG_S2F_VIS
from objects.CSCG._2d.forms.standard._2_form.base.reconstruct import _2dCSCG_S2F_Reconstruct



class _2Form_BASE(_2dCSCG_Standard_Form):
    """"""
    def __init_2form_base__(self):
        self._discretize_ = _2dCSCG_S2F_Discretize(self)
        self._visualize_ = _2dCSCG_S2F_VIS(self)
        self._reconstruct_ = None

    @property
    def visualize(self):
        return self._visualize_

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2

        if func_body.__class__.__name__ == '_2dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"2dCSCG 2form FUNC do not accept func _2dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"2dCSCG 2form FUNC do not accept func {func_body.__class__}")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2
        raise Exception(f"2dCSCG 2form BC do not accept func {func_body.__class__}")

    @property
    def discretize(self):
        return self._discretize_

    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = _2dCSCG_S2F_Reconstruct(self)
        return self._reconstruct_


    def ___PRIVATE_make_reconstruction_matrix_on_grid___(self, xi, eta):
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
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, eta)
        RM = dict()
        INDICES = self.mesh.elements.indices
        base = basis[0].T
        type_cache = dict()
        for i in INDICES:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark

            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in type_cache:
                    RM[i] = type_cache[typeWr2Metric]
                else:
                    det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                    RM_i_ = np.einsum('ji, j -> ji', base, det_iJ, optimize='greedy')
                    type_cache[typeWr2Metric] = RM_i_
                    RM[i] = RM_i_
            else:
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                RM[i] = np.einsum('ji, j -> ji', base, det_iJ, optimize='greedy')

        return RM


    def ___PRIVATE_operator_inner___(self, _, i, xietasigma, quad_weights, bfSelf, bfOther):
        """Note that here we only return local matrices."""

        element = self.mesh.elements[i]
        detJ = element.coordinate_transformation.Jacobian(*xietasigma)

        Mi = np.einsum('im, jm, m -> ij',
            bfOther[0], bfSelf[0], np.reciprocal(detJ)*quad_weights,
            optimize='greedy'
        )
        Mi = spspa.csr_matrix(Mi)
        return Mi

    def ___PRIVATE_operator_wedge___(self, other, quad_degree=None):
        """In fact, it is integral over wedge product."""
        assert other.k == 0, "Need a _2dCSCG_0Form"
        assert self.mesh == other.mesh, "Meshes do not match."
        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], other.dqp[i]])) for i in range(2)]
        quad_nodes, _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        _, basisS = self.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        _, basisO = other.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        W = np.einsum('im, jm, m -> ij', basisO[0], basisS[0], quad_weights, optimize='greedy')
        return spspa.csc_matrix(W)