# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from scipy import sparse as spspa
from _3dCSCG.forms.standard._3_form.discretize.main import _3dCSCG_Discretize
from _3dCSCG.forms.standard.base.main import _3dCSCG_Standard_Form
from _3dCSCG.forms.standard._3_form.special import _3Form_Special


class _3dCSCG_3Form(_3dCSCG_Standard_Form):
    """Standard 3-form.

    :param mesh:
    :param space:
    :param is_hybrid:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, is_hybrid=True,
        orientation='outer', numbering_parameters='Naive',  name=None):
        if name is None:
            if is_hybrid:
                name = 'hybrid-' + orientation + '-oriented-3-form'
            else:
                name = orientation + '-oriented-3-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 3
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_3form')
        self._special_ = _3Form_Special(self)
        self.___PRIVATE_reset_cache___()
        self._discretize_ = _3dCSCG_Discretize(self)
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 3form FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 3form FUNC do not accept func {func_body.__class__}")

    @property
    def special(self):
        return self._special_

    @property
    def discretize(self):
        return self._discretize_


    def reconstruct(self, xi, eta, sigma, ravel=False, i=None, regions=None):
        """
        Reconstruct the standard 3-form.

        Given ``xi``, ``eta`` and ``sigma``, we reconstruct the 3-form on ``meshgrid(xi, eta, sigma)``
        in all elements.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param sigma: A 1d iterable object of floats between -1 and 1.
        :param i: (`default`:``None``) Do the reconstruction for ``#i`` element. if it is ``None``,
            then do it for all elements.
        :type i: int, None
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :type sigma: list, tuple, numpy.ndarray
        :param bool ravel: (`default`:``False``) If we return 1d data?
        :param regions: Higher priority than input ``i``.
        :returns: A tuple of outputs

            1. (Dict[int, list]) -- :math:`x, y, z` coordinates.
            2. (Dict[int, list]) -- Reconstructed values.
        """
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, eta, sigma)
        xyz = dict()
        value = dict()
        shape = [len(xi), len(eta), len(sigma)]

        if regions is None:
            INDICES = self.mesh.elements.indices if i is None else [i, ]
        else:
            if regions == 'all':
                regions = self.mesh.domain.regions
            elif isinstance(regions, str):
                regions = [regions,]
            else:
                pass
            assert isinstance(regions, (list, tuple)), f"regions={regions} is wrong."
            assert len(set(regions)) == len(regions), f"regions={regions} has repeated regions."
            for i, r in enumerate(regions):
                assert r in self.mesh.domain.regions, f"regions[{i}]={r} is wrong."

            INDICES = list()
            for i in self.mesh.elements.indices:
                ri = self.mesh.do.FIND_region_name_of_element(i)
                if ri in regions:
                    INDICES.append(i)

        biJC = dict()
        for i in INDICES:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            if typeWr2Metric in biJC:
                basis_det_iJ = biJC[typeWr2Metric]
            else:
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                basis_det_iJ = basis[0] * det_iJ
                if isinstance(typeWr2Metric, str):
                    biJC[typeWr2Metric] = basis_det_iJ
            v = np.einsum('ij, i -> j', basis_det_iJ, self.cochain.local[i], optimize='greedy')
            if ravel:
                value[i] = [v,]
            else:
                # noinspection PyUnresolvedReferences
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                value[i] = [v.reshape(shape, order='F'),]
        return xyz, value

    def ___PRIVATE_make_reconstruction_matrix_on_grid___(self, xi, et, sg):
        """Make the reconstruction matrices for all mesh elements. These matrices are stored in
        a dict whose keys are the numbers of mesh elements and values are the local reconstruction
        matrices.

        Let `RM` be the reconstruction matrix (or the tuple of three matrices).
        If we want to do the local reconstruction, we do

            RM[i] @ f.cochain.local[i]

        and we will get the reconstructions of the form `f` on `meshgrid(xi, eta, sigma)` in mesh-element
        #i. And if `f` is a scalar form, we get a 1d array. And if `f` is a vector form, we get a
        tuple of three 1d arrays (its three components along x, y, z directions.)

        :param xi: 1d array
        :param et: 1d array
        :param sg: 1d array
        :return:
        """
        xietasigma, basis = self.do.evaluate_basis_at_meshgrid(xi, et, sg)
        biJC = dict()
        RM = dict()
        for i in self.mesh.elements:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in biJC:
                    RM[i] = biJC[typeWr2Metric]
                else:
                    det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                    basis_det_iJ = (basis[0] * det_iJ).T
                    biJC[typeWr2Metric] = basis_det_iJ
                    RM[i] = basis_det_iJ
            else:
                det_iJ = element.coordinate_transformation.inverse_Jacobian(*xietasigma)
                basis_det_iJ = basis[0] * det_iJ
                RM[i] = basis_det_iJ.T
        return RM



    def ___PRIVATE_operator_inner___(self, _, i, xietasigma, quad_weights, bfSelf, bfOther):
        """Note that here we only return a local matrix."""
        element = self.mesh.elements[i]
        detJ = element.coordinate_transformation.Jacobian(*xietasigma)
        Mi = np.einsum('im, jm, m -> ij',
            bfOther[0], bfSelf[0], np.reciprocal(detJ)*quad_weights,
            optimize='greedy'
        )
        Mi = spspa.csc_matrix(Mi)
        return Mi

    def ___PRIVATE_operator_wedge___(self, other, quad_degree=None):
        """In fact, it is integral over wedge product."""
        assert other.__class__.__name__ == '_0Form', "Need a _3dCSCG_0Form"
        assert self.mesh == other.mesh, "Meshes do not match."
        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], other.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        _, basisS = self.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        _, basisO = other.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        W = np.einsum('im, jm, m -> ij', basisO[0], basisS[0], quad_weights, optimize='greedy')
        return spspa.csc_matrix(W)






if __name__ == '__main__':
    # mpiexec python _3dCSCG\forms\standard\_3_form\main.py

    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',2), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    f3 = FC('3-f', is_hybrid=False)

    f3.TW.func.do.set_func_body_as(es, 'pressure')
    f3.TW.current_time = 0
    f3.TW.___DO_push_all_to_instant___()
    f3.do.discretize()

    # from tools.CSCG.partial_dofs import PartialDofs
    #
    # pd = PartialDofs(f3)