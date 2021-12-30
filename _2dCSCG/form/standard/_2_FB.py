
from root.config import *
from scipy import sparse as spspa
from _2dCSCG.form.standard.main import _2dCSCG_Standard_Form
from SCREWS.quadrature import Quadrature




class _2Form_BASE(_2dCSCG_Standard_Form):
    """"""
    def RESET_cache(self):
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        super().RESET_cache()

    def ___TW_FUNC_body_checker___(self, func_body):
        assert func_body.__class__.__name__ == '_2dCSCG_ScalarField'
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2


    def discretize(self, update_cochain=True, **kwargs):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param kwargs: Keyword arguments to be passed to the particular discretize method.
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if self.func.ftype == 'standard':
            return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain, **kwargs)
        else:
            raise NotImplementedError()

    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain:bool=True, quad_degree=None):
        p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None \
            or quad_degree != self.___DISCRETIZE_STANDARD_CACHE___[0]:
            magic_factor = 0.25
            xi = np.zeros((self.NUM_basis, p[0]+1, p[1]+1))
            et = np.zeros((self.NUM_basis, p[0]+1, p[1]+1))
            volume = np.zeros(self.NUM_basis)
            for j in range(self.p[1]):
                for i in range(self.p[0]):
                    m = i + j*self.p[0]
                    xi[m,...] = (quad_nodes[0][:,np.newaxis].repeat(p[1]+1, axis=1) + 1)\
                              * (self.space.nodes[0][i+1]-self.space.nodes[0][i]
                              )/2 + self.space.nodes[0][i]
                    et[m,...] = (quad_nodes[1][np.newaxis,:].repeat(p[0]+1, axis=0) + 1)\
                              * (self.space.nodes[1][j+1]-self.space.nodes[1][j]
                              )/2 + self.space.nodes[1][j]
                    volume[m] = (self.space.nodes[0][i+1]-self.space.nodes[0][i]) \
                              * (self.space.nodes[1][j+1]-self.space.nodes[1][j])  * magic_factor
            self.___DISCRETIZE_STANDARD_CACHE___ = quad_degree, xi, et, volume
        else:
            xi, et, volume = self.___DISCRETIZE_STANDARD_CACHE___[1:]

        cochainLocal = dict()
        f = self.func.body[0]
        JC = dict()
        for i in self.mesh.elements.indices:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            xyz = element.coordinate_transformation.mapping(xi, et)
            if typeWr2Metric in JC:
                detJ = JC[typeWr2Metric]
            else:
                detJ = element.coordinate_transformation.Jacobian(xi, et)
                if isinstance(typeWr2Metric, str):
                    JC[typeWr2Metric] = detJ
            fxyz = f(*xyz)
            cochainLocal[i] = np.einsum('jkl, k, l, j -> j',
                fxyz*detJ, quad_weights[0], quad_weights[1],
                volume, optimize='greedy')

        # isKronecker? ...
        if not self.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: self.cochain.local = cochainLocal
        # ...
        return cochainLocal

    def reconstruct(self, xi, eta, ravel=False, i=None):
        """
        Reconstruct the standard 3-form.

        Given ``xi``, ``eta`` and ``sigma``, we reconstruct the 3-form on ``meshgrid(xi, eta, sigma)``
        in all elements.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param i: (`default`:``None``) Do the reconstruction for ``#i`` element. if it is ``None``,
            then do it for all elements.
        :type i: int, None
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :param bool ravel: (`default`:``False``) If we return 1d data?
        :returns: A tuple of outputs

            1. (Dict[int, list]) -- :math:`x, y, z` coordinates.
            2. (Dict[int, list]) -- Reconstructed values.
        """
        xietasigma, basis = self.DO.evaluate_basis_at_meshgrid(xi, eta)
        xyz = dict()
        value = dict()
        shape = [len(xi), len(eta)]
        INDICES = self.mesh.elements.indices if i is None else [i, ]
        iJC = dict()
        for i in INDICES:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            if typeWr2Metric in iJC:
                det_iJ = iJC[typeWr2Metric]
            else:
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                if isinstance(typeWr2Metric, str):
                    iJC[typeWr2Metric] = det_iJ
            v = np.einsum('ij, i -> j', basis[0], self.cochain.local[i], optimize='greedy') * det_iJ
            if ravel:
                value[i] = [v,]
            else:
                # noinspection PyUnresolvedReferences
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                value[i] = [v.reshape(shape, order='F'),]
        return xyz, value



    def ___OPERATORS_inner___(self, _, i, xietasigma, quad_weights, bfSelf, bfOther):
        """Note that here we only return a local matrix."""

        element = self.mesh.elements[i]
        detJ = element.coordinate_transformation.Jacobian(*xietasigma)

        Mi = np.einsum('im, jm, m -> ij',
            bfOther[0], bfSelf[0], np.reciprocal(detJ)*quad_weights,
            optimize='greedy'
        )
        Mi = spspa.csr_matrix(Mi)
        return Mi

    def ___OPERATORS_wedge___(self, other, quad_degree=None):
        """In fact, it is integral over wedge product."""
        assert other.k == 0, "Need a _2dCSCG_0Form"
        assert self.mesh == other.mesh, "Meshes do not match."
        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], other.dqp[i]])) for i in range(2)]
        quad_nodes, _, quad_weights = self.space.DO_evaluate_quadrature(quad_degree)
        _, basisS = self.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        _, basisO = other.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        W = np.einsum('im, jm, m -> ij', basisO[0], basisS[0], quad_weights, optimize='greedy')
        return spspa.csc_matrix(W)