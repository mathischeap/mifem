# -*- coding: utf-8 -*-
import numpy as np
from objects.CSCG._2d.forms.standard.base.main import _2dCSCG_Standard_Form
from components.quadrature import Quadrature
from objects.CSCG._2d.forms.standard._1_form.base.visualize.main import _2dCSCG_S1F_VIS


class _1Form_BASE(_2dCSCG_Standard_Form):
    """"""
    def __init_1form_base__(self):
        self._visualize_ = _2dCSCG_S1F_VIS(self)

    @property
    def visualize(self):
        return self._visualize_


    def ___PRIVATE_do_evaluate_basis_at_meshgrid___(self, xi, eta, compute_xietasigma=True):
        """
        Evaluate the basis functions on ``meshgrid(xi, eta, sigma)``.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :param bool compute_xietasigma: (`default`:``True``) If we compute the
            ``meshgrid(xi, eta, sigma, indexing='ij')``.
        :returns: A tuple of outputs:

            1. (None, tuple) -- ``(xi, eta, sigma)`` after ``meshgrid`` and ``ravel('F')``.
            2. tuple -- The evaluated basis functions.
        """
        return self.space.do.evaluate_form_basis_at_meshgrid(
            self.k, xi, eta, orientation=self.orientation, compute_xietasigma=compute_xietasigma)


    def RESET_cache(self):
        super().RESET_cache()

    def ___Pr_check_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2

        if func_body.__class__.__name__ == '_2dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"2dCSCG 1form FUNC do not accept func _2dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"2dCSCG 1form FUNC do not accept func {func_body.__class__}")


    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2

        if func_body.__class__.__name__ == '_2dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"2dCSCG 1form BC do not accept func _2dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"2dCSCG 1form BC do not accept func {func_body.__class__}")

    def ___PRIVATE_discretize_preparation___(self, d_='', quad_degree=None):
        p = [self.dqp[i] + 3 for i in range(self.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
        p_x, p_y = p
        edges_size = [self.space.nodes[i][1:] - self.space.nodes[i][:-1] for i in range(2)]
        cell_nodes = [(0.5 * (edges_size[i][np.newaxis, :]) * (quad_nodes[i][:, np.newaxis] + 1)
            + self.space.nodes[i][:-1]).ravel('F') for i in range(2)]

        if d_ == 'x':
            quad_xi = np.tile(cell_nodes[0], self.p[1] + 1).reshape(
                (p_x + 1, self.p[0] * (self.p[1] + 1)), order='F')
            quad_eta = np.repeat(self.space.nodes[1][np.newaxis, :], self.p[0], axis=0).ravel('F')
            quad_eta = quad_eta[np.newaxis, :].repeat(p_x + 1, axis=0)
            ES = np.tile(edges_size[0], self.p[1] + 1)
            return quad_xi, quad_eta, ES, quad_weights
        elif d_ == 'y':
            quad_xi = np.tile(self.space.nodes[0], self.p[1])[np.newaxis, :].repeat(p_y + 1, axis=0)
            quad_eta = np.repeat(cell_nodes[1].reshape(
                (p_y + 1, self.p[1]), order='F'), self.p[0] + 1, axis=1)
            ES = np.repeat(edges_size[1], self.p[0] + 1)
            return quad_xi, quad_eta, ES, quad_weights
        else:
            raise Exception()