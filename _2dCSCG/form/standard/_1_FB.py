
import numpy as np
from _2dCSCG.form.standard.main import _2dCSCG_Standard_Form
from SCREWS.quadrature import Quadrature


class _1Form_BASE(_2dCSCG_Standard_Form):
    """"""

    def ___DO_evaluate_basis_at_meshgrid___(self, xi, eta, compute_xietasigma=True):
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
        return self.space.DO_evaluate_form_basis_at_meshgrid(
            self.k, xi, eta, orientation=self.orientation, compute_xietasigma=compute_xietasigma)


    def RESET_cache(self):
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        super().RESET_cache()

    def ___TW_FUNC_body_checker___(self, func_body):
        assert func_body.__class__.__name__ == '_2dCSCG_VectorField'
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2

    def discretize(self, update_cochain=True, **kwargs):
        """
        Discretize the current function (a vector field:
        :class:`_3dCSCG.form.continuous.vector._3dCSCG_VectorField`) to cochain.
        It is actually a wrapper of multiple methods that discretize functions of different types (a vector
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if self.func.ftype == 'standard':
            # noinspection PyUnresolvedReferences
            return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain, **kwargs)
        else:
            raise NotImplementedError()

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