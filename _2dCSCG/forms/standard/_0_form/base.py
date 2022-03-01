
from root.config import *
from scipy import sparse as spspa
from _2dCSCG.forms.standard.base.main import _2dCSCG_Standard_Form




class _0Form_BASE(_2dCSCG_Standard_Form):
    """"""
    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.__class__.__name__ == '_2dCSCG_ScalarField'
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def discretize(self, update_cochain=True):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if self.func.ftype == 'standard':
            return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain)
        else:
            raise NotImplementedError()

    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain=True):
        nodes = list(np.meshgrid(*self.space.nodes, indexing='ij'))
        nodes = [nodes[i].ravel('F') for i in range(2)]
        cochainLocal = dict()
        for i in self.mesh.elements:
            element = self.mesh.elements[i]
            xyz = element.coordinate_transformation.mapping(*nodes)
            cochainLocal[i] = self.func.body[0](*xyz)
        # isKronecker? ...
        if not self.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: self.cochain.local = cochainLocal
        # ...
        return cochainLocal

    def reconstruct(self, xi, eta, ravel=False, i=None):
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, eta)
        xyz = dict()
        value = dict()
        shape = [len(xi), len(eta)]
        INDICES = self.mesh.elements.indices if i is None else [i, ]
        for i in INDICES:
            element = self.mesh.elements[i]
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            v = np.einsum('ij, i -> j', basis[0], self.cochain.local[i], optimize='optimal')
            if ravel:
                value[i] = [v,]
            else:
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                value[i] = [v.reshape(shape, order='F'),]
        return xyz, value

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
        quad_nodes , _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        xietasigma, basisS = self.do.evaluate_basis_at_meshgrid(*quad_nodes)
        _, basisO = other.do.evaluate_basis_at_meshgrid(*quad_nodes)
        W = np.einsum('im, jm, m -> ij', basisO[0], basisS[0], quad_weights, optimize='greedy')
        return spspa.csc_matrix(W)