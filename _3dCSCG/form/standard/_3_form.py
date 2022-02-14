# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config import *
from scipy import sparse as spspa
from SCREWS.frozen import FrozenOnly
from SCREWS.quadrature import Quadrature
from _3dCSCG.form.standard.main import _3dCSCG_Standard_Form


class _3Form(_3dCSCG_Standard_Form):
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
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        super().RESET_cache()

    def ___TW_FUNC_body_checker___(self, func_body):
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
            if self.TW.func.body.__class__.__name__ == '_3dCSCG_ScalarField':

                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain, **kwargs)
                else:
                    raise NotImplementedError(f"3dCSCG 3-form cannot (target func) discretize _3dCSCG_ScalarField of ftype={self.func.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 3-form can not (target func) discretize {self.TW.func.body.__class__}.')
        else:
            raise NotImplementedError(f"3dCSCG 3-form cannot discretize while targeting at {target}.")

    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain:bool=True, quad_degree=None):
        p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None \
            or quad_degree != self.___DISCRETIZE_STANDARD_CACHE___[0]:
            magic_factor = 0.125
            xi = np.zeros((self.NUM_basis, p[0]+1, p[1]+1, p[2]+1))
            et = np.zeros((self.NUM_basis, p[0]+1, p[1]+1, p[2]+1))
            si = np.zeros((self.NUM_basis, p[0]+1, p[1]+1, p[2]+1))
            volume = np.zeros(self.NUM_basis)
            for k in range(self.p[2]):
                for j in range(self.p[1]):
                    for i in range(self.p[0]):
                        m = i + j*self.p[0] + k*self.p[0]*self.p[1]
                        xi[m,...] = (quad_nodes[0][:,np.newaxis].repeat(p[1]+1,
                          axis=1)[:,:,np.newaxis].repeat(p[2]+1, axis=2) + 1)\
                                  * (self.space.nodes[0][i+1]-self.space.nodes[0][i]
                                  )/2 + self.space.nodes[0][i]
                        et[m,...] = (quad_nodes[1][np.newaxis,:].repeat(p[0]+1,
                          axis=0)[:,:,np.newaxis].repeat(p[2]+1, axis=2) + 1)\
                                  * (self.space.nodes[1][j+1]-self.space.nodes[1][j]
                                  )/2 + self.space.nodes[1][j]
                        si[m,...] = (quad_nodes[2][np.newaxis,:].repeat(p[1]+1,
                          axis=0)[np.newaxis,:,:].repeat(p[0]+1, axis=0) + 1)\
                                  * (self.space.nodes[2][k+1]-self.space.nodes[2][k]
                                  )/2 + self.space.nodes[2][k]
                        volume[m] = (self.space.nodes[0][i+1]-self.space.nodes[0][i]) \
                                  * (self.space.nodes[1][j+1]-self.space.nodes[1][j]) \
                                  * (self.space.nodes[2][k+1]-self.space.nodes[2][k]) * magic_factor
            self.___DISCRETIZE_STANDARD_CACHE___ = quad_degree, xi, et, si, volume
        else:
            xi, et, si, volume = self.___DISCRETIZE_STANDARD_CACHE___[1:]
        cochainLocal = dict()
        f = self.func.body[0]
        JC = dict()
        for i in self.mesh.elements.indices:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            xyz = element.coordinate_transformation.mapping(xi, et, si)
            if typeWr2Metric in JC:
                detJ = JC[typeWr2Metric]
            else:
                detJ = element.coordinate_transformation.Jacobian(xi, et, si)
                if isinstance(typeWr2Metric, str):
                    JC[typeWr2Metric] = detJ
            fxyz = f(*xyz)
            cochainLocal[i] = np.einsum('jklm, k, l, m, j -> j',
                fxyz*detJ, quad_weights[0], quad_weights[1], quad_weights[2],
                volume, optimize='greedy'
            )
        # isKronecker? ...
        if not self.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: self.cochain.local = cochainLocal
        # ...
        return cochainLocal

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
        xietasigma, basis = self.DO.evaluate_basis_at_meshgrid(xi, eta, sigma)
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
                ri = self.mesh.DO.FIND_region_name_of_element(i)
                if ri in regions:
                    INDICES.append(i)

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
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
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
        Mi = spspa.csc_matrix(Mi)
        return Mi

    def ___OPERATORS_wedge___(self, other, quad_degree=None):
        """In fact, it is integral over wedge product."""
        assert other.__class__.__name__ == '_0Form', "Need a _3dCSCG_0Form"
        assert self.mesh == other.mesh, "Meshes do not match."
        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], other.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = self.space.DO_evaluate_quadrature(quad_degree)
        _, basisS = self.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        _, basisO = other.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        W = np.einsum('im, jm, m -> ij', basisO[0], basisS[0], quad_weights, optimize='greedy')
        return spspa.csc_matrix(W)


class _3Form_Special(FrozenOnly):
    def __init__(self, _3sf):
        self._sf_ = _3sf
        self._freeze_self_()




if __name__ == '__main__':
    # mpiexec python _3dCSCG\form\standard\_3_form.py

    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',2), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    f3 = FC('3-f', is_hybrid=False)

    f3.TW.func.DO.set_func_body_as(es, 'pressure')
    f3.TW.current_time = 0
    f3.TW.___DO_push_all_to_instant___()
    f3.DO.discretize()

    from TOOLS.CSCG.partial_dofs import PartialDofs

    pd = PartialDofs(f3)