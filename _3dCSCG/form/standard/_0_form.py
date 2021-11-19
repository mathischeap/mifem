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
from _3dCSCG.form.standard.main import _3dCSCG_Standard_Form


class _0Form(_3dCSCG_Standard_Form):
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
        if name is None: name = orientation + '-oriented-0-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_0form')
        self._special_ = _0Form_Special(self)
        self.RESET_cache()
        self._freeze_self_()

    def ___TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0form FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form FUNC do not accept func {func_body.__class__}")

    def ___TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 0form BC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form BC do not accept func {func_body.__class__}")


    def RESET_cache(self):
        super().RESET_cache()

    @property
    def special(self):
        return self._special_

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
            if self.TW.func.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain)
                else:
                    raise NotImplementedError(f"3dCSCG 0-form cannot (target func) discretize _3dCSCG_ScalarField of ftype={self.func.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 0-form can not (target func) discretize {self.TW.func.body.__class__}.')

        elif target == 'BC':
            if self.TW.BC.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if self.BC.ftype == 'standard':
                    # always do not update cochain & and target always be "BC"
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=False, target='BC')

                elif self.BC.ftype == "boundary-wise":
                    # we will always not update cochain & and always set target to be "BC"
                    return self.___PRIVATE_discretize_boundary_wise_ftype___()

                else:
                    raise NotImplementedError(f"3dCSCG 0-form cannot (target BC) discretize _3dCSCG_ScalarField of ftype={self.BC.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 0-form can not (target BC) discretize {self.TW.BC.body.__class__}.')
        else:
            raise NotImplementedError(f"3dCSCG 0-form cannot discretize while targeting at {target}.")


    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain=True, target='func'):
        nodes = list(np.meshgrid(*self.space.nodes, indexing='ij'))
        nodes = [nodes[i].ravel('F') for i in range(3)]
        cochainLocal = dict()
        if target == 'func':
            FUNC = self.func.body[0]
        elif target == 'BC':
            FUNC = self.BC.body[0]
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"_0Form.___PRIVATE_discretize_standard_ftype___ "
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

    def ___PRIVATE_discretize_boundary_wise_ftype___(self):
        """

        :return:
        """
        nodes = list(np.meshgrid(*self.space.nodes, indexing='ij'))
        nodes = [nodes[i].ravel('F') for i in range(3)]
        FUNC = self.BC.body
        RANGE_element_sides = self.mesh.boundaries.RANGE_element_sides
        cochainLocal = dict()
        for bn in FUNC:
            func_bn = FUNC[bn][0]
            element_sides = RANGE_element_sides[bn]
            elements, sides = list(), list()
            for element_side in element_sides:
                element = int(element_side[:-1])
                side = element_side[-1]
                elements.append(element)
                sides.append(side)
            for i, side in zip(elements, sides):
                element = self.mesh.elements[i]
                xyz = element.coordinate_transformation.mapping(*nodes)
                local_cochain = func_bn(*xyz)
                local_dofs = self.numbering.DO.FIND.local_dofs_on_element_side(side)
                if i not in cochainLocal:
                    cochainLocal[i] = dict()

                cochainLocal[i][side] = local_cochain[local_dofs]

        return 'Boundary only local cochain', cochainLocal





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

        for i in INDICES:
            element = self.mesh.elements[i]
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            v = np.einsum('ij, i -> j', basis[0], self.cochain.local[i], optimize='optimal')
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
        Mi = np.einsum('im, jm, m -> ij', bfOther[0], bfSelf[0], detJ*quad_weights, optimize='greedy')
        Mi = spspa.csc_matrix(Mi)
        return Mi




class _0Form_Special(FrozenOnly):
    def __init__(self, _0sf):
        self._sf_ = _0sf
        self._freeze_self_()





if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\standard\_0_form.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([5,5,5])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    # es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    def p(t, x, y, z):
        return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t


    scalar = FC('scalar', p)

    f0 = FC('0-f', is_hybrid=True)
    f0.TW.func.DO.set_func_body_as(scalar)
    f0.TW.current_time = 0
    f0.TW.DO.push_all_to_instant()
    f0.DO.discretize()
    print(f0.error.L())


    df0 = FC('0-adf')
    df0.prime.TW.func.DO.set_func_body_as(scalar)
    df0.prime.TW.current_time = 0
    df0.prime.TW.DO.push_all_to_instant()
    df0.prime.DO.discretize()
    print(df0.prime.error.L())


    # region = mesh.domain.regions['R:R']
    # RS = region.sub_geometry.GENERATE_perpendicular_slice_object(t=4.5/9)
    # MS = mesh.sub_geometry.GENERATE_perpendicular_slice_object(RS)
    # f0.visualize.matplot.perpendicular_slice(MS, saveto='')

    # f0.export.field.to_file('t0f.mat')
    # print(f0.matrices.trace)