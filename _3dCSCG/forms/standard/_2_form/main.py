# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('../')

from root.config import *
from scipy import sparse as spspa
from screws.quadrature import Quadrature
from _3dCSCG.forms.standard.base.main import _3dCSCG_Standard_Form
from _3dCSCG.forms.standard._2_form.special.main import _2Form_Special
from _3dCSCG.forms.standard._2_form.project.main import _2Form_Projection


class _3dCSCG_2Form(_3dCSCG_Standard_Form):
    """
    Standard 2-form.

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
                name = 'hybrid-' + orientation + '-oriented-2-form'
            else:
                name = orientation + '-oriented-2-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_2form')
        self._special_ = _2Form_Special(self)
        self._projection_ = _2Form_Projection(self)
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 2form FUNC do not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 2form FUNC do not accept func {func_body.__class__}")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 2form BC do not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 2form BC do not accept func {func_body.__class__}")


    @property
    def special(self):
        return self._special_

    @property
    def projection(self):
        """A wrapper of all projection methods."""
        return self._projection_

    def discretize(self, update_cochain=True, target='func', **kwargs):
        """
        Discretize the current function (a vector field:
        :class:`_3dCSCG.form.continuous.vector._3dCSCG_VectorField`) to cochain.
        It is actually a wrapper of multiple methods that discretize functions of different types (a vector
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if target == 'func':

            if self.TW.func.body.__class__.__name__ == '_3dCSCG_VectorField':

                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain, **kwargs)

                else:
                    raise NotImplementedError(f"3dCSCG 2-form cannot (target func) discretize _3dCSCG_VectorField of ftype={self.func.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 2-form can not (target func) discretize {self.TW.func.body.__class__}.')

        elif target == 'BC':
            if self.TW.BC.body.__class__.__name__ == '_3dCSCG_VectorField':

                if self.BC.ftype == 'standard':
                    # always do not update cochain & and target always be "BC"
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=False, target='BC', **kwargs)

                elif self.BC.ftype == "boundary-wise":
                    # we will always not update cochain & and always set target to be "BC"
                    return self.___PRIVATE_discretize_boundary_wise_ftype___(**kwargs)

                else:
                    raise NotImplementedError(f"3dCSCG 2-form cannot (target BC) discretize _3dCSCG_VectorField of ftype={self.BC.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 2-form can not (target BC) discretize {self.TW.BC.body.__class__}.')

        else:
            raise NotImplementedError(f"3dCSCG 2-form cannot discretize while targeting at {target}.")

    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain=True, target='func', quad_degree=None):
        """The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.

        :param update_cochain:
        :param target:
        :param quad_degree:
        :return:
        """
        p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None or \
                quad_degree != self.___DISCRETIZE_STANDARD_CACHE___['quadDegree']:
            self.___DISCRETIZE_STANDARD_CACHE___ = dict()

            xi = np.zeros((self.NUM_basis_components[0], p[1] + 1, p[2] + 1))
            et = np.zeros((self.NUM_basis_components[0], p[1] + 1, p[2] + 1))
            si = np.zeros((self.NUM_basis_components[0], p[1] + 1, p[2] + 1))
            area_dydz = np.zeros((self.NUM_basis_components[0]))
            for k in range(self.p[2]):
                for j in range(self.p[1]):
                    for i in range(self.p[0] + 1):
                        m = i + j * (self.p[0] + 1) + k * (self.p[0] + 1) * self.p[1]
                        xi[m, ...] = np.ones((p[1] + 1, p[2] + 1)) * self.space.nodes[0][i]
                        et[m, ...] = (quad_nodes[1][:, np.newaxis].repeat(p[2] + 1, axis=1) + 1) * (
                                self.space.nodes[1][j + 1] - self.space.nodes[1][j]) / 2 + \
                                     self.space.nodes[1][j]
                        si[m, ...] = (quad_nodes[2][np.newaxis, :].repeat((p[1] + 1), axis=0) + 1) * (
                                self.space.nodes[2][k + 1] - self.space.nodes[2][k]) / 2 + \
                                     self.space.nodes[2][k]
                        area_dydz[m] = (self.space.nodes[2][k + 1] - self.space.nodes[2][k]) \
                                       * (self.space.nodes[1][j + 1] - self.space.nodes[1][j])

            self.___DISCRETIZE_STANDARD_CACHE___['X'] = (xi, et, si, area_dydz)

            xi = np.zeros((self.NUM_basis_components[1], p[0] + 1, p[2] + 1))
            et = np.zeros((self.NUM_basis_components[1], p[0] + 1, p[2] + 1))
            si = np.zeros((self.NUM_basis_components[1], p[0] + 1, p[2] + 1))
            area_dzdx = np.zeros((self.NUM_basis_components[1]))
            for k in range(self.p[2]):
                for j in range(self.p[1] + 1):
                    for i in range(self.p[0]):
                        m = i + j * self.p[0] + k * (self.p[1] + 1) * self.p[0]
                        xi[m, ...] = (quad_nodes[0][:, np.newaxis].repeat(p[2] + 1, axis=1) + 1) * (
                                self.space.nodes[0][i + 1] - self.space.nodes[0][i]) / 2 + \
                                     self.space.nodes[0][i]
                        et[m, ...] = np.ones((p[0] + 1, p[2] + 1)) * self.space.nodes[1][j]
                        si[m, ...] = (quad_nodes[2][np.newaxis, :].repeat(p[0] + 1, axis=0) + 1) * (
                                self.space.nodes[2][k + 1] - self.space.nodes[2][k]) / 2 + \
                                     self.space.nodes[2][k]
                        area_dzdx[m] = (self.space.nodes[2][k + 1] - self.space.nodes[2][k]) \
                                       * (self.space.nodes[0][i + 1] - self.space.nodes[0][i])

            self.___DISCRETIZE_STANDARD_CACHE___['Y'] = (xi, et, si, area_dzdx)

            xi = np.zeros((self.NUM_basis_components[2], p[0] + 1, p[1] + 1))
            et = np.zeros((self.NUM_basis_components[2], p[0] + 1, p[1] + 1))
            si = np.zeros((self.NUM_basis_components[2], p[0] + 1, p[1] + 1))
            area_dxdy = np.zeros((self.NUM_basis_components[2]))
            for k in range(self.p[2] + 1):
                for j in range(self.p[1]):
                    for i in range(self.p[0]):
                        m = i + j * self.p[0] + k * self.p[1] * self.p[0]
                        xi[m, ...] = (quad_nodes[0][:, np.newaxis].repeat(p[1] + 1, axis=1) + 1) * (
                                self.space.nodes[0][i + 1] - self.space.nodes[0][i]) / 2 + \
                                     self.space.nodes[0][i]
                        et[m, ...] = (quad_nodes[1][np.newaxis, :].repeat(p[0] + 1, axis=0) + 1) * (
                                self.space.nodes[1][j + 1] - self.space.nodes[1][j]) / 2 + \
                                     self.space.nodes[1][j]
                        si[m, ...] = np.ones((p[0] + 1, p[1] + 1)) * self.space.nodes[2][k]
                        area_dxdy[m] = (self.space.nodes[1][j + 1] - self.space.nodes[1][j]) \
                                       * (self.space.nodes[0][i + 1] - self.space.nodes[0][i])

            self.___DISCRETIZE_STANDARD_CACHE___['Z'] = (xi, et, si, area_dxdy)
            self.___DISCRETIZE_STANDARD_CACHE___['quadDegree'] = quad_degree
            del xi, et, si
        else:
            pass

        xi_x, et_x, si_x, area_dydz = self.___DISCRETIZE_STANDARD_CACHE___['X']
        xi_y, et_y, si_y, area_dzdx = self.___DISCRETIZE_STANDARD_CACHE___['Y']
        xi_z, et_z, si_z, area_dxdy = self.___DISCRETIZE_STANDARD_CACHE___['Z']

        local_dydz = dict()
        local_dzdx = dict()
        local_dxdy = dict()

        JXC, JYC, JZC = dict(), dict(), dict()

        if target == 'func':
            FUNC = self.func
        elif target == 'BC':
            FUNC = self.BC
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise NotImplementedError(
                f"_2Form.___PRIVATE_discretize_standard_ftype___ "
                f"does not work for target={target}.")

        for i in self.mesh.elements:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark

            smctm_x = element.coordinate_transformation.mapping(xi_x, et_x, si_x)
            if typeWr2Metric in JXC:
                Jx_0, Jx_1, Jx_2 = JXC[typeWr2Metric]
            else:
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    J11 = element.coordinate_transformation.J11(xi_x, et_x, si_x)
                    J22 = element.coordinate_transformation.J22(xi_x, et_x, si_x)
                    Jx_0 = J11 * J22
                    Jx_1 = 0
                    Jx_2 = 0
                else:
                    J = element.coordinate_transformation.Jacobian_matrix(xi_x, et_x, si_x)
                    Jx_0 = J[1][1] * J[2][2] - J[1][2] * J[2][1]
                    Jx_1 = J[2][1] * J[0][2] - J[2][2] * J[0][1]
                    Jx_2 = J[0][1] * J[1][2] - J[0][2] * J[1][1]
                if isinstance(typeWr2Metric, str):
                    JXC[typeWr2Metric] = Jx_0, Jx_1, Jx_2

            smctm_y = element.coordinate_transformation.mapping(xi_y, et_y, si_y)
            if typeWr2Metric in JYC:
                Jy_0, Jy_1, Jy_2 = JYC[typeWr2Metric]
            else:
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    J00 = element.coordinate_transformation.J00(xi_y, et_y, si_y)
                    J22 = element.coordinate_transformation.J22(xi_y, et_y, si_y)
                    Jy_0 = 0
                    Jy_1 = J00 * J22
                    Jy_2 = 0
                else:
                    J = element.coordinate_transformation.Jacobian_matrix(xi_y, et_y, si_y)
                    Jy_0 = J[1][2] * J[2][0] - J[1][0] * J[2][2]
                    Jy_1 = J[2][2] * J[0][0] - J[2][0] * J[0][2]
                    Jy_2 = J[0][2] * J[1][0] - J[0][0] * J[1][2]
                if isinstance(typeWr2Metric, str):
                    JYC[typeWr2Metric] = Jy_0, Jy_1, Jy_2

            smctm_z = element.coordinate_transformation.mapping(xi_z, et_z, si_z)
            if typeWr2Metric in JZC:
                Jz_0, Jz_1, Jz_2 = JZC[typeWr2Metric]
            else:
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    J11 = element.coordinate_transformation.J11(xi_z, et_z, si_z)
                    J00 = element.coordinate_transformation.J00(xi_z, et_z, si_z)
                    Jz_0 = 0
                    Jz_1 = 0
                    Jz_2 = J11 * J00
                else:
                    J = element.coordinate_transformation.Jacobian_matrix(xi_z, et_z, si_z)
                    Jz_0 = J[1][0] * J[2][1] - J[1][1] * J[2][0]
                    Jz_1 = J[2][0] * J[0][1] - J[2][1] * J[0][0]
                    Jz_2 = J[0][0] * J[1][1] - J[0][1] * J[1][0]
                if isinstance(typeWr2Metric, str):
                    JZC[typeWr2Metric] = Jz_0, Jz_1, Jz_2

            # Now it is time to do the reduction.

            u = FUNC.body[0](*smctm_x)
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                uvw_dydz = Jx_0 * u
            else:
                v = FUNC.body[1](*smctm_x)
                w = FUNC.body[2](*smctm_x)
                uvw_dydz = Jx_0 * u + Jx_1 * v + Jx_2 * w
            local_dydz[i] = self.___PRIVATE_discretize_standard_einsum___(
                uvw_dydz, quad_weights[1], quad_weights[2], area_dydz
            )
            v = FUNC.body[1](*smctm_y)
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                uvw_dzdx = Jy_1 * v
            else:
                u = FUNC.body[0](*smctm_y)
                w = FUNC.body[2](*smctm_y)
                uvw_dzdx = Jy_0 * u + Jy_1 * v + Jy_2 * w
            local_dzdx[i] = self.___PRIVATE_discretize_standard_einsum___(
                uvw_dzdx, quad_weights[0], quad_weights[2], area_dzdx
            )
            w = FUNC.body[2](*smctm_z)
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                uvw_dxdy = Jz_2 * w
            else:
                u = FUNC.body[0](*smctm_z)
                v = FUNC.body[1](*smctm_z)
                uvw_dxdy = Jz_0 * u + Jz_1 * v + Jz_2 * w
            local_dxdy[i] = self.___PRIVATE_discretize_standard_einsum___(
                uvw_dxdy, quad_weights[0], quad_weights[1], area_dxdy
            )
        del JXC, JYC, JZC
        # isKronecker? ...
        if not self.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        cochainLocal = dict()
        for i in self.mesh.elements.indices:
            cochainLocal[i] = np.hstack((local_dydz[i], local_dzdx[i], local_dxdy[i]))
        if update_cochain: self.cochain.local = cochainLocal
        # 'locally full local cochain': provide cochain.local and locally they are full for all local dofs
        return 'locally full local cochain', cochainLocal


    def ___PRIVATE_discretize_boundary_wise_ftype___(self, quad_degree=None):
        """
        'Boundary only local cochain' means we return a dict, its keys are mesh-element numbers,
        its values are also dictionaries whose keys are mesh-element-side names, like 'N', 'S' and
        so on, and values are the mesh-element-side(trace-element)-wise local cochains. For example
        cochainLocal = {
                1: {'N': [4, 3, 1, 1.5, ...], 'W': [...]},
                23: {...},
                ...}
        We know we have cochains for mesh-element #1, #23, ..., and for mesh element #1, we have
        local cochain on its North side and West side.

        :param quad_degree:
        :return:
        """
        p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None or \
                quad_degree != self.___DISCRETIZE_STANDARD_CACHE___['quadDegree']:
            self.___DISCRETIZE_STANDARD_CACHE___ = dict()

            xi = np.zeros((self.NUM_basis_components[0], p[1] + 1, p[2] + 1))
            et = np.zeros((self.NUM_basis_components[0], p[1] + 1, p[2] + 1))
            si = np.zeros((self.NUM_basis_components[0], p[1] + 1, p[2] + 1))
            area_dydz = np.zeros((self.NUM_basis_components[0]))
            for k in range(self.p[2]):
                for j in range(self.p[1]):
                    for i in range(self.p[0] + 1):
                        m = i + j * (self.p[0] + 1) + k * (self.p[0] + 1) * self.p[1]
                        xi[m, ...] = np.ones((p[1] + 1, p[2] + 1)) * self.space.nodes[0][i]
                        et[m, ...] = (quad_nodes[1][:, np.newaxis].repeat(p[2] + 1, axis=1) + 1) * (
                                self.space.nodes[1][j + 1] - self.space.nodes[1][j]) / 2 + \
                                     self.space.nodes[1][j]
                        si[m, ...] = (quad_nodes[2][np.newaxis, :].repeat((p[1] + 1), axis=0) + 1) * (
                                self.space.nodes[2][k + 1] - self.space.nodes[2][k]) / 2 + \
                                     self.space.nodes[2][k]
                        area_dydz[m] = (self.space.nodes[2][k + 1] - self.space.nodes[2][k]) \
                                       * (self.space.nodes[1][j + 1] - self.space.nodes[1][j])

            self.___DISCRETIZE_STANDARD_CACHE___['X'] = (xi, et, si, area_dydz)

            xi = np.zeros((self.NUM_basis_components[1], p[0] + 1, p[2] + 1))
            et = np.zeros((self.NUM_basis_components[1], p[0] + 1, p[2] + 1))
            si = np.zeros((self.NUM_basis_components[1], p[0] + 1, p[2] + 1))
            area_dzdx = np.zeros((self.NUM_basis_components[1]))
            for k in range(self.p[2]):
                for j in range(self.p[1] + 1):
                    for i in range(self.p[0]):
                        m = i + j * self.p[0] + k * (self.p[1] + 1) * self.p[0]
                        xi[m, ...] = (quad_nodes[0][:, np.newaxis].repeat(p[2] + 1, axis=1) + 1) * (
                                self.space.nodes[0][i + 1] - self.space.nodes[0][i]) / 2 + \
                                     self.space.nodes[0][i]
                        et[m, ...] = np.ones((p[0] + 1, p[2] + 1)) * self.space.nodes[1][j]
                        si[m, ...] = (quad_nodes[2][np.newaxis, :].repeat(p[0] + 1, axis=0) + 1) * (
                                self.space.nodes[2][k + 1] - self.space.nodes[2][k]) / 2 + \
                                     self.space.nodes[2][k]
                        area_dzdx[m] = (self.space.nodes[2][k + 1] - self.space.nodes[2][k]) \
                                       * (self.space.nodes[0][i + 1] - self.space.nodes[0][i])

            self.___DISCRETIZE_STANDARD_CACHE___['Y'] = (xi, et, si, area_dzdx)

            xi = np.zeros((self.NUM_basis_components[2], p[0] + 1, p[1] + 1))
            et = np.zeros((self.NUM_basis_components[2], p[0] + 1, p[1] + 1))
            si = np.zeros((self.NUM_basis_components[2], p[0] + 1, p[1] + 1))
            area_dxdy = np.zeros((self.NUM_basis_components[2]))
            for k in range(self.p[2] + 1):
                for j in range(self.p[1]):
                    for i in range(self.p[0]):
                        m = i + j * self.p[0] + k * self.p[1] * self.p[0]
                        xi[m, ...] = (quad_nodes[0][:, np.newaxis].repeat(p[1] + 1, axis=1) + 1) * (
                                self.space.nodes[0][i + 1] - self.space.nodes[0][i]) / 2 + \
                                     self.space.nodes[0][i]
                        et[m, ...] = (quad_nodes[1][np.newaxis, :].repeat(p[0] + 1, axis=0) + 1) * (
                                self.space.nodes[1][j + 1] - self.space.nodes[1][j]) / 2 + \
                                     self.space.nodes[1][j]
                        si[m, ...] = np.ones((p[0] + 1, p[1] + 1)) * self.space.nodes[2][k]
                        area_dxdy[m] = (self.space.nodes[1][j + 1] - self.space.nodes[1][j]) \
                                       * (self.space.nodes[0][i + 1] - self.space.nodes[0][i])

            self.___DISCRETIZE_STANDARD_CACHE___['Z'] = (xi, et, si, area_dxdy)
            self.___DISCRETIZE_STANDARD_CACHE___['quadDegree'] = quad_degree
            del xi, et, si
        else:
            pass

        xi_x, et_x, si_x, area_dydz = self.___DISCRETIZE_STANDARD_CACHE___['X']
        xi_y, et_y, si_y, area_dzdx = self.___DISCRETIZE_STANDARD_CACHE___['Y']
        xi_z, et_z, si_z, area_dxdy = self.___DISCRETIZE_STANDARD_CACHE___['Z']

        FUNC = self.BC.body

        JXC, JYC, JZC = dict(), dict(), dict()
        RANGE_element_sides = self.mesh.boundaries.range_of_element_sides

        cochainLocal = dict()
        for bn in FUNC:
            func_bn = FUNC[bn]
            element_sides = RANGE_element_sides[bn]
            elements, sides = list(), list()
            for element_side in element_sides:
                element = int(element_side[:-1])
                side = element_side[-1]
                elements.append(element)
                sides.append(side)

            for i, side in zip(elements, sides):
                element = self.mesh.elements[i]
                typeWr2Metric = element.type_wrt_metric.mark

                smctm_x = element.coordinate_transformation.mapping(xi_x, et_x, si_x)
                if typeWr2Metric in JXC:
                    Jx_0, Jx_1, Jx_2 = JXC[typeWr2Metric]
                else:
                    if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                        J11 = element.coordinate_transformation.J11(xi_x, et_x, si_x)
                        J22 = element.coordinate_transformation.J22(xi_x, et_x, si_x)
                        Jx_0 = J11 * J22
                        Jx_1 = 0
                        Jx_2 = 0
                    else:
                        J = element.coordinate_transformation.Jacobian_matrix(xi_x, et_x, si_x)
                        Jx_0 = J[1][1] * J[2][2] - J[1][2] * J[2][1]
                        Jx_1 = J[2][1] * J[0][2] - J[2][2] * J[0][1]
                        Jx_2 = J[0][1] * J[1][2] - J[0][2] * J[1][1]
                    if isinstance(typeWr2Metric, str):
                        JXC[typeWr2Metric] = Jx_0, Jx_1, Jx_2

                smctm_y = element.coordinate_transformation.mapping(xi_y, et_y, si_y)
                if typeWr2Metric in JYC:
                    Jy_0, Jy_1, Jy_2 = JYC[typeWr2Metric]
                else:
                    if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                        J00 = element.coordinate_transformation.J00(xi_y, et_y, si_y)
                        J22 = element.coordinate_transformation.J22(xi_y, et_y, si_y)
                        Jy_0 = 0
                        Jy_1 = J00 * J22
                        Jy_2 = 0
                    else:
                        J = element.coordinate_transformation.Jacobian_matrix(xi_y, et_y, si_y)
                        Jy_0 = J[1][2] * J[2][0] - J[1][0] * J[2][2]
                        Jy_1 = J[2][2] * J[0][0] - J[2][0] * J[0][2]
                        Jy_2 = J[0][2] * J[1][0] - J[0][0] * J[1][2]
                    if isinstance(typeWr2Metric, str):
                        JYC[typeWr2Metric] = Jy_0, Jy_1, Jy_2

                smctm_z = element.coordinate_transformation.mapping(xi_z, et_z, si_z)
                if typeWr2Metric in JZC:
                    Jz_0, Jz_1, Jz_2 = JZC[typeWr2Metric]
                else:
                    if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                        J11 = element.coordinate_transformation.J11(xi_z, et_z, si_z)
                        J00 = element.coordinate_transformation.J00(xi_z, et_z, si_z)
                        Jz_0 = 0
                        Jz_1 = 0
                        Jz_2 = J11 * J00
                    else:
                        J = element.coordinate_transformation.Jacobian_matrix(xi_z, et_z, si_z)
                        Jz_0 = J[1][0] * J[2][1] - J[1][1] * J[2][0]
                        Jz_1 = J[2][0] * J[0][1] - J[2][1] * J[0][0]
                        Jz_2 = J[0][0] * J[1][1] - J[0][1] * J[1][0]
                    if isinstance(typeWr2Metric, str):
                        JZC[typeWr2Metric] = Jz_0, Jz_1, Jz_2

                # Now it is time to do the reduction.....................
                u = func_bn[0](*smctm_x)
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    uvw_dydz = Jx_0 * u
                else:
                    v = func_bn[1](*smctm_x)
                    w = func_bn[2](*smctm_x)
                    uvw_dydz = Jx_0 * u + Jx_1 * v + Jx_2 * w
                LOCAL_FULL_COCHAIN_x = self.___PRIVATE_discretize_standard_einsum___(
                    uvw_dydz, quad_weights[1], quad_weights[2], area_dydz
                )
                v = func_bn[1](*smctm_y)
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    uvw_dzdx = Jy_1 * v
                else:
                    u = func_bn[0](*smctm_y)
                    w = func_bn[2](*smctm_y)
                    uvw_dzdx = Jy_0 * u + Jy_1 * v + Jy_2 * w
                LOCAL_FULL_COCHAIN_y = self.___PRIVATE_discretize_standard_einsum___(
                    uvw_dzdx, quad_weights[0], quad_weights[2], area_dzdx
                )
                w = func_bn[2](*smctm_z)
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    uvw_dxdy = Jz_2 * w
                else:
                    u = func_bn[0](*smctm_z)
                    v = func_bn[1](*smctm_z)
                    uvw_dxdy = Jz_0 * u + Jz_1 * v + Jz_2 * w
                LOCAL_FULL_COCHAIN_z = self.___PRIVATE_discretize_standard_einsum___(
                    uvw_dxdy, quad_weights[0], quad_weights[1], area_dxdy
                )

                LOCAL_FULL_COCHAIN = np.hstack((LOCAL_FULL_COCHAIN_x,
                                                LOCAL_FULL_COCHAIN_y,
                                                LOCAL_FULL_COCHAIN_z))

                local_dofs = self.numbering.do.find.local_dofs_on_element_side(side)
                if i not in cochainLocal:
                    cochainLocal[i] = dict()

                cochainLocal[i][side] = LOCAL_FULL_COCHAIN[local_dofs]
        # 'Boundary only local cochain': provide cochain.local and only for local dofs on the element sides.
        return 'Boundary only local cochain', cochainLocal





    @staticmethod
    def ___PRIVATE_discretize_standard_einsum___(uvw, quad_weights_1, quad_weights_2, area):
        """ """
        return np.einsum('jkl, kl, j -> j',
                         uvw, np.tensordot(quad_weights_1, quad_weights_2, axes=0),
                         area * 0.25, optimize='optimal'
                         )

    def reconstruct(self, xi, eta, sigma, ravel=False, i=None, regions=None):
        """
        Do the reconstruction.

        :param xi:
        :param eta:
        :param sigma:
        :param ravel:
        :param i: The element we want to reconstruct in. When it is ``None``, we do the reconstruction for
            all elements and store the results in one coordinate dictionary and one value dictionary.
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
                ri = self.mesh.do.FIND_region_name_of_element(i)
                if ri in regions:
                    INDICES.append(i)

        iJ = self.mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
        JC = dict() # local cache
        for i in INDICES:
            element = self.mesh.elements[i]
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            u = np.einsum('ij, i -> j', basis[0], self.cochain.local_('x')[i], optimize='greedy')
            v = np.einsum('ij, i -> j', basis[1], self.cochain.local_('y')[i], optimize='greedy')
            w = np.einsum('ij, i -> j', basis[2], self.cochain.local_('z')[i], optimize='greedy')
            value[i] = [None, None, None]
            typeWr2Metric = element.type_wrt_metric.mark

            if typeWr2Metric in JC:
                _0u, _0v, _0w, _1u, _1v, _1w, _2u, _2v, _2w = JC[typeWr2Metric]
            else:
                iJi = iJ[i]
                if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                    _0u = iJi[1][1] * iJi[2][2]
                    _0v = None
                    _0w = None
                    _1u = None
                    _1v = iJi[2][2] * iJi[0][0]
                    _1w = None
                    _2u = None
                    _2v = None
                    _2w = iJi[0][0] * iJi[1][1]
                else:
                    _0u = iJi[1][1] * iJi[2][2] - iJi[1][2] * iJi[2][1]
                    _0v = iJi[2][1] * iJi[0][2] - iJi[2][2] * iJi[0][1]
                    _0w = iJi[0][1] * iJi[1][2] - iJi[0][2] * iJi[1][1]
                    _1u = iJi[1][2] * iJi[2][0] - iJi[1][0] * iJi[2][2]
                    _1v = iJi[2][2] * iJi[0][0] - iJi[2][0] * iJi[0][2]
                    _1w = iJi[0][2] * iJi[1][0] - iJi[0][0] * iJi[1][2]
                    _2u = iJi[1][0] * iJi[2][1] - iJi[1][1] * iJi[2][0]
                    _2v = iJi[2][0] * iJi[0][1] - iJi[2][1] * iJi[0][0]
                    _2w = iJi[0][0] * iJi[1][1] - iJi[0][1] * iJi[1][0]

                if isinstance(typeWr2Metric, str):
                    JC[typeWr2Metric] = _0u, _0v, _0w, _1u, _1v, _1w, _2u, _2v, _2w

            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                value[i][0] = u * _0u
                value[i][1] = v * _1v
                value[i][2] = w * _2w
            else:
                value[i][0] = u * _0u +  v * _0v +  w * _0w
                value[i][1] = u * _1u +  v * _1v +  w * _1w
                value[i][2] = u * _2u +  v * _2v +  w * _2w

            if ravel:
                pass
            else:
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                # noinspection PyUnresolvedReferences
                value[i] = [value[i][j].reshape(shape, order='F') for j in range(3)]
        return xyz, value


    def ___PRIVATE_operator_inner___(self, other, i, xietasigma, quad_weights, bfSelf, bfOther):
        """
        We compute the inner product between ``self`` and ``other`` in element ``i``.

        Note that here we only return a local matrix.
        """
        element = self.mesh.elements[i]
        mark = element.type_wrt_metric.mark
        J = element.coordinate_transformation.Jacobian_matrix(*xietasigma)
        sqrtg = element.coordinate_transformation.Jacobian(*xietasigma, J=J)
        iJ = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=J)
        g = element.coordinate_transformation.inverse_metric_matrix(*xietasigma, iJ=iJ)
        del J, iJ
        if isinstance(mark, str) and mark[:4] == 'Orth':
            M00 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[1][1]*g[2][2], bfOther[0], bfSelf[0])
            M11 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[2][2]*g[0][0], bfOther[1], bfSelf[1])
            M22 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[0][0]*g[1][1], bfOther[2], bfSelf[2])
            M01 = None
            M02 = None
            M12 = None
            M10 = None
            M20 = None
            M21 = None
        else:
            M00 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[1][1]*g[2][2]-g[1][2]*g[2][1], bfOther[0], bfSelf[0])
            M11 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[2][2]*g[0][0]-g[2][0]*g[0][2], bfOther[1], bfSelf[1])
            M22 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[0][0]*g[1][1]-g[0][1]*g[1][0], bfOther[2], bfSelf[2])
            g12_20_g10_22 = g[1][2] * g[2][0] - g[1][0] * g[2][2]
            g10_21_g11_20 = g[1][0] * g[2][1] - g[1][1] * g[2][0]
            g20_01_g21_00 = g[2][0] * g[0][1] - g[2][1] * g[0][0]
            M01 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g12_20_g10_22, bfOther[0], bfSelf[1])
            M02 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g10_21_g11_20, bfOther[0], bfSelf[2])
            M12 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g20_01_g21_00, bfOther[1], bfSelf[2])
            if other is self:
                M10 = M01.T
                M20 = M02.T
                M21 = M12.T
            else:
                M10 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g12_20_g10_22, bfOther[1], bfSelf[0])
                M20 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g10_21_g11_20, bfOther[2], bfSelf[0])
                M21 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g20_01_g21_00, bfOther[2], bfSelf[1])
        Mi = spspa.bmat([(M00, M01, M02),
                         (M10, M11, M12),
                         (M20, M21, M22)], format='csc')
        return Mi

    @staticmethod
    def ___PRIVATE_inner_H1___(quad_weights, sqrtg, g, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights*sqrtg*g, bfO, bfS, optimize='optimal')
        return spspa.csc_matrix(M)






if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\form\standard\_2_form.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.0)([14,14,14])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)


    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    velocity = FC('vector', (u,v,w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)


    f2 = FC('2-f', is_hybrid=False)


    f2.TW.func.do.set_func_body_as(velocity)
    f2.TW.current_time = 0
    f2.TW.___DO_push_all_to_instant___()
    f2.discretize()

    SO = f2.special.vortex_detection.Q_and_lambda2([-1, 0, 0.5, 1], [-1, 1], [-1, -0.5, 0, 0.5, 1])

    #
    # f0x, f0y, f0z = f2.projection.to.vector_of_3_standard_0forms()
    #
    #
    # f0x.TW.func.do.set_func_body_as(U)
    # f0x.TW.current_time = 0
    # f0x.TW.do.push_all_to_instant()
    # print(f0x.error.L())
    # f0y.TW.func.do.set_func_body_as(V)
    # f0y.TW.current_time = 0
    # f0y.TW.do.push_all_to_instant()
    # print(f0y.error.L())
    # f0z.TW.func.do.set_func_body_as(W)
    # f0z.TW.current_time = 0
    # f0z.TW.do.push_all_to_instant()
    # print(f0z.error.L())
