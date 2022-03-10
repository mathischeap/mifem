
from root.config.main import *
from scipy import sparse as spspa
from _2dCSCG.forms.standard.base.main import _2dCSCG_Standard_Form
from screws.quadrature import Quadrature




class _2Form_BASE(_2dCSCG_Standard_Form):
    """"""
    def ___PRIVATE_reset_cache___(self):
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2

        if func_body.__class__.__name__ == '_2dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"2dCSCG 2form FUNC do not accept func _2dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 2form FUNC do not accept func {func_body.__class__}")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3
        raise Exception(f"3dCSCG 2form BC do not accept func {func_body.__class__}")


    def discretize(self, update_cochain=True, target='func', **kwargs):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :param kwargs: Keyword arguments to be passed to the particular discretize method.
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if target == 'func':
            if self.func.ftype == 'standard':
                return self.___PRIVATE_discretize_standard_ftype___(
                    update_cochain=update_cochain, target=target,**kwargs)
            else:
                raise NotImplementedError()
        elif target == 'BC':
            raise NotImplementedError(f'2dCSCG 1-form can not (target BC) '
                                      f'discretize {self.TW.BC.body.__class__}.')
        else:
            raise NotImplementedError(f"2dCSCG 1-form cannot discretize "
                                      f"while targeting at {target}.")


    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain:bool=True, target='func', quad_degree=None):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.

        :param update_cochain:
        :param quad_degree:
        :return:
        """
        p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None \
            or quad_degree != self.___DISCRETIZE_STANDARD_CACHE___[0]:
            magic_factor = 0.25
            xi = np.zeros((self.num.basis, p[0]+1, p[1]+1))
            et = np.zeros((self.num.basis, p[0]+1, p[1]+1))
            volume = np.zeros(self.num.basis)
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

        # ----------- target ---------------------------------------------------
        cochainLocal = dict()
        if target == 'func':
            FUNC = self.func.body[0]
        else:
            raise NotImplementedError(f"I cannot deal with target = {target}.")
        # ============================================================================

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
            fxyz = FUNC(*xyz)
            cochainLocal[i] = np.einsum('jkl, k, l, j -> j',
                fxyz*detJ, quad_weights[0], quad_weights[1],
                volume, optimize='greedy')

        # isKronecker? ...
        if not self.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: self.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal

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
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, eta)
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
                basis_det_iJ = iJC[typeWr2Metric]
            else:
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                basis_det_iJ = basis[0] * det_iJ
                if isinstance(typeWr2Metric, str):
                    iJC[typeWr2Metric] = basis_det_iJ

            v = np.einsum('ij, i -> j', basis_det_iJ, self.cochain.local[i], optimize='greedy')
            if ravel:
                value[i] = [v,]
            else:
                # noinspection PyUnresolvedReferences
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(2)]
                value[i] = [v.reshape(shape, order='F'),]
        return xyz, value

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