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
from TOOLS.linear_algebra.elementwise_cache import EWC_SparseMatrix
from _3dCSCG.form.standard._0_form import _0Form


class _2Form(_3dCSCG_Standard_Form):
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
        if name is None: name = orientation + '-oriented-2-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_2form')
        self._special_ = _2Form_Special(self)
        self._projection_ = _2Form_Projection(self)
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        super().RESET_cache()

    def ___TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 2form FUNC do not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 2form FUNC do not accept func {func_body.__class__}")

    def ___TW_BC_body_checker___(self, func_body):
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
        """"""
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
        RANGE_element_sides = self.mesh.boundaries.RANGE_element_sides

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

                local_dofs = self.numbering.DO.FIND.local_dofs_on_element_side(side)
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


    def ___OPERATORS_inner___(self, other, i, xietasigma, quad_weights, bfSelf, bfOther):
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




class _2Form_Special(FrozenOnly):
    def __init__(self, _2sf):
        self._sf_ = _2sf
        self._vortex_detection_ = None
        self._freeze_self_()

    def cross_product(self, u, e, quad_degree=None):
        """
        We do ``(self X other, e)`` where ``self`` and ``other`` both are n-form, n be either 1 or 2.

        :return:
        """
        SCP_generator = ___3dCSCG_2Form_CrossProduct___(self._sf_, u, e, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')

    @property
    def vortex_detection(self):
        if self._vortex_detection_ is None:
            self._vortex_detection_ = ___3dCSCG_2Form_Vortex_Detection___(self._sf_)
        return self._vortex_detection_


class ___3dCSCG_2Form_CrossProduct___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(a \\times b, e)`.

    :param a:
    :param b:
    :param e:
    :param quad_degree:
    """
    def __init__(self, a, b, e, quad_degree=None):
        assert a.ndim == b.ndim == e.ndim, " <_3dCSCG_2StddFm_InnerWedgeMatrix> "
        assert a.k == b.k == e.k == 2, " <_3dCSCG_2StddFm_InnerWedgeMatrix> "
        assert a.mesh == b.mesh, "Meshes do not match."
        assert a.mesh == e.mesh, "Meshes do not match."
        self._mesh_ = a.mesh
        self._a_ = a
        self._b_ = b
        self._e_ = e

        if quad_degree is None:
            quad_degree = [int(np.max([a.dqp[i], b.dqp[i], e.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = a.space.DO_evaluate_quadrature(quad_degree)
        xietasigma, abf = a.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=True)
        _         , bbf = b.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        _         , ebf = e.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        self._qw_ = quad_weights
        self._abf_ = abf
        self._bbf_ = bbf
        self._ebf_ = ebf

        self._JM_    = self._mesh_.elements.coordinate_transformation.Jacobian_matrix(*xietasigma)
        self._sqrtg_ = self._mesh_.elements.coordinate_transformation.Jacobian(*xietasigma, J=self._JM_)
        self._iJ_    = self._mesh_.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=self._JM_)
        self._g_     = self._mesh_.elements.coordinate_transformation.inverse_metric_matrix(*xietasigma, iJ=self._iJ_)

        self.DO_reset_cache()
        self._freeze_self_()

    def DO_reset_cache(self):
        self._J_cache_ = dict()
        self._G_cache_ = dict()

    def _J_(self, i):
        element = self._mesh_.elements[i]
        typeWr2Metric = element.type_wrt_metric.mark
        if typeWr2Metric in self._J_cache_:
            return self._J_cache_[typeWr2Metric]
        else:
            iJ = self._iJ_[i]
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                J00 = iJ[1][1] * iJ[2][2]
                J01 = 0
                J02 = 0
                J10 = 0
                J11 = iJ[2][2] * iJ[0][0]
                J12 = 0
                J20 = 0
                J21 = 0
                J22 = iJ[0][0] * iJ[1][1]
            else:
                J00 = iJ[1][1] * iJ[2][2] - iJ[1][2] * iJ[2][1]
                J01 = iJ[2][1] * iJ[0][2] - iJ[2][2] * iJ[0][1]
                J02 = iJ[0][1] * iJ[1][2] - iJ[0][2] * iJ[1][1]
                J10 = iJ[1][2] * iJ[2][0] - iJ[1][0] * iJ[2][2]
                J11 = iJ[2][2] * iJ[0][0] - iJ[2][0] * iJ[0][2]
                J12 = iJ[0][2] * iJ[1][0] - iJ[0][0] * iJ[1][2]
                J20 = iJ[1][0] * iJ[2][1] - iJ[1][1] * iJ[2][0]
                J21 = iJ[2][0] * iJ[0][1] - iJ[2][1] * iJ[0][0]
                J22 = iJ[0][0] * iJ[1][1] - iJ[0][1] * iJ[1][0]
            J = (J00, J01, J02, J10, J11, J12, J20, J21, J22)
            # cache it even for unique mesh cells (because we may use them multiple times when do temporal iterations.)
            self._J_cache_[typeWr2Metric] = J
            return J

    def _G_(self, i):
        element = self._mesh_.elements[i]
        typeWr2Metric = element.type_wrt_metric.mark
        if typeWr2Metric in self._G_cache_:
            return self._G_cache_[typeWr2Metric]
        else:
            sqrtg = self._sqrtg_[i]
            g = self._g_[i]
            JM = self._JM_[i]
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                G00 = sqrtg * g[1][1] * g[2][2] * JM[1][1] * JM[2][2]
                G01 = 0
                G02 = 0
                G10 = 0
                G11 = sqrtg * g[2][2] * g[0][0] * JM[2][2] * JM[0][0]
                G12 = 0
                G20 = 0
                G21 = 0
                G22 = sqrtg * g[0][0] * g[1][1] * JM[0][0] * JM[1][1]
            else:
                G00 = sqrtg * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
                G01 = sqrtg * (g[1][2] * g[2][0] - g[1][0] * g[2][2])
                G02 = sqrtg * (g[1][0] * g[2][1] - g[1][1] * g[2][0])
                G10 = G01
                G11 = sqrtg * (g[2][2] * g[0][0] - g[2][0] * g[0][2])
                G12 = sqrtg * (g[2][0] * g[0][1] - g[2][1] * g[0][0])
                G20 = G02
                G21 = G12
                G22 = sqrtg * (g[0][0] * g[1][1] - g[0][1] * g[1][0])
            G = (G00, G01, G02, G10, G11, G12, G20, G21, G22)
            # cache it even for unique mesh cells (because we may use them multiple times when do temporal iterations.)
            self._G_cache_[typeWr2Metric] = G
            return G

    def __call__(self, i):
        typeWr2Metric = self._mesh_.elements[i].type_wrt_metric.mark

        a0, a1, a2 =  self._abf_ # a; given
        b0, b1, b2 =  self._bbf_ # b
        e0, e1, e2 =  self._ebf_ # epsilon

        a0p = np.einsum('ij, i -> j', a0, self._a_.cochain.local_('x')[i], optimize='greedy')
        a1p = np.einsum('ij, i -> j', a1, self._a_.cochain.local_('y')[i], optimize='greedy')
        a2p = np.einsum('ij, i -> j', a2, self._a_.cochain.local_('z')[i], optimize='greedy')

        J00, J01, J02, J10, J11, J12, J20, J21, J22 = self._J_(i)
        G00, G01, G02, G10, G11, G12, G20, G21, G22 = self._G_(i)

        if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
            a0 = a0p * J00
            a1 = a1p * J11
            a2 = a2p * J22

            A01, A02 =  J00*a2, -J00*a1
            A10, A12 = -J11*a2,  J11*a0
            A20, A21 =  J22*a1, -J22*a0

            m01 = A10*G00
            m02 = A20*G00
            m10 = A01*G11
            m12 = A21*G11
            m20 = A02*G22
            m21 = A12*G22

            M01 = np.einsum('iw, jw, w -> ij', e0, b1, m01*self._qw_, optimize='greedy')
            M02 = np.einsum('iw, jw, w -> ij', e0, b2, m02*self._qw_, optimize='greedy')
            M10 = np.einsum('iw, jw, w -> ij', e1, b0, m10*self._qw_, optimize='greedy')
            M12 = np.einsum('iw, jw, w -> ij', e1, b2, m12*self._qw_, optimize='greedy')
            M20 = np.einsum('iw, jw, w -> ij', e2, b0, m20*self._qw_, optimize='greedy')
            M21 = np.einsum('iw, jw, w -> ij', e2, b1, m21*self._qw_, optimize='greedy')

            M = ([None                 , spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
                 [spspa.csc_matrix(M10), None                 , spspa.csc_matrix(M12)],
                 [spspa.csc_matrix(M20), spspa.csc_matrix(M21), None                 ])

        else:
            #TODO: correct the following ...
            raise NotImplementedError()

            # a0 = a0p * J00 + a1p * J01 + a2p * J02
            # a1 = a0p * J10 + a1p * J11 + a2p * J12
            # a2 = a0p * J20 + a1p * J21 + a2p * J22
            #
            # A00, A01, A02 = (J20*a1 - J10*a2), (J00*a2 - J20*a0), (J10*a0 - J00*a1)
            # A10, A11, A12 = (J21*a1 - J11*a2), (J01*a2 - J21*a0), (J11*a0 - J01*a1)
            # A20, A21, A22 = (J22*a1 - J12*a2), (J02*a2 - J22*a0), (J12*a0 - J02*a1)
            #
            # m00 = A00*G00 + A01*G01 + A02*G02
            # m01 = A10*G00 + A11*G01 + A12*G02
            # m02 = A20*G00 + A21*G01 + A22*G02
            # m10 = A00*G10 + A01*G11 + A02*G12
            # m11 = A10*G10 + A11*G11 + A12*G12
            # m12 = A20*G10 + A21*G11 + A22*G12
            # m20 = A00*G20 + A01*G21 + A02*G22
            # m21 = A10*G20 + A11*G21 + A12*G22
            # m22 = A20*G20 + A21*G21 + A22*G22
            #
            # M00 = np.einsum('iw, jw, w -> ij', e0, b0, m00*self._qw_, optimize='greedy')
            # M01 = np.einsum('iw, jw, w -> ij', e0, b1, m01*self._qw_, optimize='greedy')
            # M02 = np.einsum('iw, jw, w -> ij', e0, b2, m02*self._qw_, optimize='greedy')
            # M10 = np.einsum('iw, jw, w -> ij', e1, b0, m10*self._qw_, optimize='greedy')
            # M11 = np.einsum('iw, jw, w -> ij', e1, b1, m11*self._qw_, optimize='greedy')
            # M12 = np.einsum('iw, jw, w -> ij', e1, b2, m12*self._qw_, optimize='greedy')
            # M20 = np.einsum('iw, jw, w -> ij', e2, b0, m20*self._qw_, optimize='greedy')
            # M21 = np.einsum('iw, jw, w -> ij', e2, b1, m21*self._qw_, optimize='greedy')
            # M22 = np.einsum('iw, jw, w -> ij', e2, b2, m22*self._qw_, optimize='greedy')
            # M = ([spspa.csc_matrix(M00), spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
            #      [spspa.csc_matrix(M10), spspa.csc_matrix(M11), spspa.csc_matrix(M12)],
            #      [spspa.csc_matrix(M20), spspa.csc_matrix(M21), spspa.csc_matrix(M22)])

        MW = spspa.bmat(M, format='csc')
        return MW

class ___3dCSCG_2Form_Vortex_Detection___(FrozenOnly):
    """A wrapper of all vortex detection methods. So, we consider this 1 form as
    a variable of a flow field."""
    def __init__(self, _2sf):
        self._sf_ = _2sf
        self._freeze_self_()

    def ___PRIVATE_generate_gradient_tensor_at___(self, xi, eta, sigma):
        """We compute the gradient tensor of this 2form. To do so, we first project
        this 2-form into a vector of 3 standard 0-forms which represent the three
        components. Then we do the gradient (apply the incidence matrix E10) to each
        standard 0-form.

        It returns a 3 by 3 tensor representing
            ((du_dx, du_dy, du_dz),
             (dv_dx, dv_dy, dv_dz),
             (dw_dx, dw_dy, dw_dz)).
        Each value are 3d evaluated at *meshgrid(xi, eta, sigma, indexing='ij)

        :param xi: 1d increasing array in [-1,1]
        :param eta: 1d increasing array in [-1,1]
        :param sigma: 1d increasing array in [-1,1]

        """
        assert np.ndim(xi) == 1 and np.all(np.diff(xi) >0) and np.max(xi) <= 1 and np.min(xi) >= -1, \
            f"xi={xi} wrong, should be 1d array in [-1,1] and increasing."
        assert np.ndim(eta) == 1 and np.all(np.diff(eta) >0) and np.max(eta) <= 1 and np.min(eta) >= -1, \
            f"eta={eta} wrong, should be 1d array in [-1,1] and increasing."
        assert np.ndim(sigma) == 1 and np.all(np.diff(sigma) >0) and np.max(sigma) <= 1 and np.min(sigma) >= -1, \
            f"sigma={sigma} wrong, should be 1d array in [-1,1] and increasing."

        U, V, W = self._sf_.projection.to.vector_of_3_standard_0forms()
        dU = U.coboundary()
        dV = V.coboundary()
        dW = W.coboundary()

        xyz, rdU = dU.reconstruct(xi, eta, sigma, ravel=False)
        xyz, rdV = dV.reconstruct(xi, eta, sigma, ravel=False)
        xyz, rdW = dW.reconstruct(xi, eta, sigma, ravel=False)

        dU_dx, dU_dy, dU_dz = dict(), dict(), dict()
        dV_dx, dV_dy, dV_dz = dict(), dict(), dict()
        dW_dx, dW_dy, dW_dz = dict(), dict(), dict()

        for i in rdU:
            dU_dx[i], dU_dy[i], dU_dz[i] = rdU[i]
            dV_dx[i], dV_dy[i], dV_dz[i] = rdV[i]
            dW_dx[i], dW_dy[i], dW_dz[i] = rdW[i]

        return xyz, ((dU_dx, dU_dy, dU_dz),
                     (dV_dx, dV_dy, dV_dz),
                     (dW_dx, dW_dy, dW_dz))

    def ___PRIVATE_generate_S_and_Omega___(self, xi, eta, sigma):
        """
        S and Omega are the symmetric and antisymmetric components of gradient
        tensor G. So both S and Omega are 3 by 3 tensor, and

        S_{i,j} = 0.5 * (G_{i,j} + G_{j,i})
        Omega_{i,j} = 0.5 * (G_{i,j} - G_{j,i})

        Each value are 3d evaluated at *meshgrid(xi, eta, sigma, indexing='ij)

        :param xi: 1d increasing array in [-1,1]
        :param eta: 1d increasing array in [-1,1]
        :param sigma: 1d increasing array in [-1,1]

        """
        S_00, S_01, S_02 = dict(), dict(), dict()
        S_10, S_11, S_12 = dict(), dict(), dict()
        S_20, S_21, S_22 = dict(), dict(), dict()

        O_00, O_01, O_02 = dict(), dict(), dict()
        O_10, O_11, O_12 = dict(), dict(), dict()
        O_20, O_21, O_22 = dict(), dict(), dict()

        xyz, GT = self.___PRIVATE_generate_gradient_tensor_at___(xi, eta, sigma)

        dU_xyz, dV_xyz, dW_xyz = GT
        U_00, U_01, U_02 = dU_xyz
        U_10, U_11, U_12 = dV_xyz
        U_20, U_21, U_22 = dW_xyz

        for i in U_00: # will go through all local mesh elements
            u_00 = U_00[i]
            u_01 = U_01[i]
            u_02 = U_02[i]
            u_10 = U_10[i]
            u_11 = U_11[i]
            u_12 = U_12[i]
            u_20 = U_20[i]
            u_21 = U_21[i]
            u_22 = U_22[i]

            s_00 = u_00
            s_11 = u_11
            s_22 = u_22

            s_01 = 0.5 * (u_01 + u_10)
            s_02 = 0.5 * (u_02 + u_20)
            s_10 = s_01
            s_12 = 0.5 * (u_12 + u_21)
            s_20 = s_02
            s_21 = s_12
            S_00[i], S_01[i], S_02[i] = s_00, s_01, s_02
            S_10[i], S_11[i], S_12[i] = s_10, s_11, s_12
            S_20[i], S_21[i], S_22[i] = s_20, s_21, s_22

            o_00 = o_11 = o_22 = 0
            o_01 = 0.5 * (u_01 - u_10)
            o_02 = 0.5 * (u_02 - u_20)
            o_10 = - o_01
            o_12 = 0.5 * (u_12 - u_21)
            o_20 = - o_02
            o_21 = - o_12
            O_00[i], O_01[i], O_02[i] = o_00, o_01, o_02
            O_10[i], O_11[i], O_12[i] = o_10, o_11, o_12
            O_20[i], O_21[i], O_22[i] = o_20, o_21, o_22

        return xyz, ((S_00, S_01, S_02), (S_10, S_11, S_12), (S_20, S_21, S_22)), \
                    ((O_00, O_01, O_02), (O_10, O_11, O_12), (O_20, O_21, O_22))

    def ___PRIVATE_generate_lambda_1_2_3_Q___(self, xi, eta, sigma):
        """  See [on the identification of a vortex] by Jeong and Hussain.

        Compute the lambda_2, and Q definitions.

        Each value are 3d evaluated at *meshgrid(xi, eta, sigma, indexing='ij)

        :param xi: 1d increasing array in [-1,1]
        :param eta: 1d increasing array in [-1,1]
        :param sigma: 1d increasing array in [-1,1]
        """
        xyz, S, O = self.___PRIVATE_generate_S_and_Omega___(xi, eta, sigma)

        S0, S1, S2 = S
        O0, O1, O2 = O

        S00, S01, S02 = S0
        S10, S11, S12 = S1
        S20, S21, S22 = S2

        O00, O01, O02 = O0
        O10, O11, O12 = O1
        O20, O21, O22 = O2

        Q = dict()
        LAMBDA_2 = dict()

        for i in S00: # we go through all local mesh elements
            s00, s01, s02 = S00[i], S01[i], S02[i]
            s10, s11, s12 = S10[i], S11[i], S12[i]
            s20, s21, s22 = S20[i], S21[i], S22[i]
            o00 ,o01, o02 = O00[i], O01[i], O02[i]
            o10, o11, o12 = O10[i], O11[i], O12[i]
            o20, o21, o22 = O20[i], O21[i], O22[i]

            SHAPE = s00.shape

            so_00 = (s00**2 + o00**2).ravel('F')
            so_01 = (s01**2 + o01**2).ravel('F')
            so_02 = (s02**2 + o02**2).ravel('F')

            so_10 = (s10**2 + o10**2).ravel('F')
            so_11 = (s11**2 + o11**2).ravel('F')
            so_12 = (s12**2 + o12**2).ravel('F')

            so_20 = (s20**2 + o20**2).ravel('F')
            so_21 = (s21**2 + o21**2).ravel('F')
            so_22 = (s22**2 + o22**2).ravel('F')

            so = np.zeros((len(so_00),3,3))
            so[:,0,0] = so_00
            so[:,0,1] = so_01
            so[:,0,2] = so_02
            so[:,1,0] = so_10
            so[:,1,1] = so_11
            so[:,1,2] = so_12
            so[:,2,0] = so_20
            so[:,2,1] = so_21
            so[:,2,2] = so_22

            eigen_values, _ = np.linalg.eig(so)
            Q[i] = (- 0.5 * np.sum(eigen_values, axis=1)).reshape(SHAPE, order='F')
            eigen_values = np.sort(eigen_values, axis=1)

            # print(eigen_values, eigen_values[:,1] )

            LAMBDA_2[i] = eigen_values[:,1].reshape(SHAPE, order='F')

        return xyz, Q, LAMBDA_2

    def Q_and_lambda2(self, xi, eta, sigma):
        """ See [on the identification of a vortex] by Jeong and Hussain.

        Each value are 3d evaluated at *meshgrid(xi, eta, sigma, indexing='ij)

        :param xi: 1d increasing array in [-1,1]
        :param eta: 1d increasing array in [-1,1]
        :param sigma: 1d increasing array in [-1,1]
        """
        xyz, Q, LAMBDA_2 = self.___PRIVATE_generate_lambda_1_2_3_Q___(xi, eta, sigma)[:3]
        return xyz, Q, LAMBDA_2




class _2Form_Projection(FrozenOnly):
    """A wrapper of all projection methods."""
    def __init__(self,_2sf):
        self._sf_ = _2sf
        self._to_ = ___3dCSCG_2Form_Project_To___(self._sf_)
        self._van_ = ___3dCSCG_2Form_Project_Van___(self._sf_)
        self._freeze_self_()

    @property
    def to(self):
        return self._to_

    @property
    def van(self):
        return self._van_

class ___3dCSCG_2Form_Project_To___(FrozenOnly):
    """A wrapper of all projection into methods"""
    def __init__(self,_2sf):
        self._sf_ = _2sf
        self._freeze_self_()

    def vector_of_3_standard_0forms(self):
        """project this 2form into a tuple of three 0forms. Each 0form stands for
        a component of the 2form (as a vector).

        Since the 0forms with the same space will be of higher degree, so we do not
        need to make a new space of higher degree. And we of course use the same mesh.
        Thus, both the mesh and space will the same as the those of this 2form."""
        space = self._sf_.space
        mesh = self._sf_.mesh

        f0_x = _0Form(mesh, space, is_hybrid=self._sf_.IS_hybrid,
            orientation=self._sf_.orientation,
            numbering_parameters=self._sf_.numbering._numbering_parameters_,
            name='Projected_x_0form_of_'+self._sf_.standard_properties.name)

        f0_y = _0Form(mesh, space, is_hybrid=self._sf_.IS_hybrid,
            orientation=self._sf_.orientation,
            numbering_parameters=self._sf_.numbering._numbering_parameters_,
            name='Projected_y_0form_of_'+self._sf_.standard_properties.name)

        f0_z = _0Form(mesh, space, is_hybrid=self._sf_.IS_hybrid,
            orientation=self._sf_.orientation,
            numbering_parameters=self._sf_.numbering._numbering_parameters_,
            name='Projected_z_0form_of_'+self._sf_.standard_properties.name)

        _, v = self._sf_.reconstruct(*space.nodes, ravel=True)
        fx_c = dict()
        fy_c = dict()
        fz_c = dict()
        for i in v: # go thorough all local mesh elements
            vx, vy, vz = v[i] # values are actually used as the local cochains of the 0forms

            fx_c[i] = vx
            fy_c[i] = vy
            fz_c[i] = vz

        f0_x.cochain.local = fx_c
        f0_y.cochain.local = fy_c
        f0_z.cochain.local = fz_c

        return f0_x, f0_y, f0_z


class ___3dCSCG_2Form_Project_Van___(FrozenOnly):
    """A wrapper of all projection into methods"""
    def __init__(self,_2sf):
        self._sf_ = _2sf
        self._freeze_self_()


if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\form\standard\_2_form.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

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


    f2.TW.func.DO.set_func_body_as(velocity)
    f2.TW.current_time = 0
    f2.TW.___DO_push_all_to_instant___()
    f2.discretize()

    SO = f2.special.vortex_detection.Q_and_lambda2([-1, 0, 0.5, 1], [-1, 1], [-1, -0.5, 0, 0.5, 1])

    #
    # f0x, f0y, f0z = f2.projection.to.vector_of_3_standard_0forms()
    #
    #
    # f0x.TW.func.DO.set_func_body_as(U)
    # f0x.TW.current_time = 0
    # f0x.TW.DO.push_all_to_instant()
    # print(f0x.error.L())
    # f0y.TW.func.DO.set_func_body_as(V)
    # f0y.TW.current_time = 0
    # f0y.TW.DO.push_all_to_instant()
    # print(f0y.error.L())
    # f0z.TW.func.DO.set_func_body_as(W)
    # f0z.TW.current_time = 0
    # f0z.TW.DO.push_all_to_instant()
    # print(f0z.error.L())
