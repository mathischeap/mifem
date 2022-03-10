
from root.config.main import *
from scipy import sparse as spspa
from _2dCSCG.forms.standard.base.main import _2dCSCG_Standard_Form




class _0Form_BASE(_2dCSCG_Standard_Form):
    """"""
    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 2

        if func_body.__class__.__name__ == '_2dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"2dCSCG 0form FUNC do not accept func _2dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form FUNC do not accept func {func_body.__class__}")


    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3
        raise Exception(f"3dCSCG 0form BC do not accept func {func_body.__class__}")



    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def discretize(self, update_cochain=True, target='func'):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if target == 'func':
            if self.TW.func.body.__class__.__name__ == '_2dCSCG_ScalarField':
                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain)
                else:
                    raise NotImplementedError(f"2dCSCG outer 0-form cannot (target func) "
                                              f"discretize _2dCSCG_ScalarField of ftype={self.func.ftype}")

            else:
                raise NotImplementedError(f'2dCSCG outer 0-form can not (target func) '
                                          f'discretize {self.TW.func.body.__class__}.')

        elif target == 'BC':
                raise NotImplementedError(f'2dCSCG outer 0-form can not (target BC) '
                                          f'discretize {self.TW.BC.body.__class__}.')
        else:
            raise NotImplementedError(f"2dCSCG outer 0-form cannot discretize "
                                      f"while targeting at {target}.")

    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain=True, target='func'):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.
        """
        nodes = list(np.meshgrid(*self.space.nodes, indexing='ij'))
        nodes = [nodes[i].ravel('F') for i in range(2)]
        cochainLocal = dict()

        if target == 'func':
            FUNC = self.func.body[0]
        elif target == 'BC':
            FUNC = self.BC.body[0]
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"_2d_outer_0Form.___PRIVATE_discretize_standard_ftype___ "
                f"does not work for target={target}.")

        for i in self.mesh.elements:
            element = self.mesh.elements[i]
            xyz = element.coordinate_transformation.mapping(*nodes)
            cochainLocal[i] = FUNC(*xyz)
        # isKronecker? ...
        if not self.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: self.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal

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
        _, basis = self.do.evaluate_basis_at_meshgrid(xi, eta)
        RM = dict()
        INDICES = self.mesh.elements.indices
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
        quad_nodes , _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        xietasigma, basisS = self.do.evaluate_basis_at_meshgrid(*quad_nodes)
        _, basisO = other.do.evaluate_basis_at_meshgrid(*quad_nodes)
        W = np.einsum('im, jm, m -> ij', basisO[0], basisS[0], quad_weights, optimize='greedy')
        return spspa.csc_matrix(W)