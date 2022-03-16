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
from _3dCSCG.forms.standard.base.main import _3dCSCG_Standard_Form
from _3dCSCG.forms.standard._0_form.special.main import _0Form_Special
from _3dCSCG.forms.standard._0_form.discretize.main import _3dCSCG_Discretize





class _3dCSCG_0Form(_3dCSCG_Standard_Form):
    """
    Standard 0-form.

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
                name = 'hybrid-' + orientation + '-oriented-0-form'
            else:
                name = orientation + '-oriented-0-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_0form')
        self._special_ = _0Form_Special(self)
        self.___PRIVATE_reset_cache___()
        self._discretize_ = _3dCSCG_Discretize(self)
        self._freeze_self_()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0form FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form FUNC do not accept func {func_body.__class__}")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 0form BC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form BC do not accept func {func_body.__class__}")


    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    @property
    def special(self):
        return self._special_

    @property
    def discretize(self):
        return self._discretize_





    def reconstruct(self, xi, eta, sigma, ravel=False, i=None, regions=None):
        """

        :param xi:
        :param eta:
        :param sigma:
        :param ravel:
        :param i:
        :param regions: Higher priority than input i.
        :return:
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
                ri = self.mesh.do.find.region_name_of_element(i)
                if ri in regions:
                    INDICES.append(i)

        for i in INDICES:
            element = self.mesh.elements[i]
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            v = np.einsum('ij, i -> j', basis[0], self.cochain.local[i], optimize='optimal')
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
        _, basis = self.do.evaluate_basis_at_meshgrid(xi, et, sg)
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









if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\forms\standard\_0_form\main.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([5,5,5])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    # es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    def p(t, x, y, z):
        return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t


    scalar = FC('scalar', p)

    f0 = FC('0-f', is_hybrid=True)
    f0.TW.func.do.set_func_body_as(scalar)
    f0.TW.current_time = 0
    f0.TW.do.push_all_to_instant()
    f0.do.discretize()
    print(f0.error.L())


    df0 = FC('0-adf')
    df0.prime.TW.func.do.set_func_body_as(scalar)
    df0.prime.TW.current_time = 0
    df0.prime.TW.do.push_all_to_instant()
    df0.prime.do.discretize()
    print(df0.prime.error.L())

    #
    # print(f0.___parameters___)

    # regions = mesh.domain.regions['R:R']
    # RS = regions.sub_geometry.make_a_perpendicular_slice_object_on(t=4.5/9)
    # MS = mesh.sub_geometry.make_a_perpendicular_slice_object_on(RS)
    # f0.visualize.matplot.perpendicular_slice(MS, saveto='')

    # f0.export.field.to_file('t0f.mat')
    # print(f0.matrices.trace)