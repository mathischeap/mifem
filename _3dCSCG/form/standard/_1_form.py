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
from TOOLS.linear_algebra.data_structures import GlobalVector
from _3dCSCG.form.standard._0_form import _0Form
from _3dCSCG.form.standard._2_form import _2Form


class _1Form(_3dCSCG_Standard_Form):
    """
    Standard 1-form.

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
                name = 'hybrid-' + orientation + '-oriented-1-form'
            else:
                name = orientation + '-oriented-1-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_1form')
        self._special_ = _1Form_Special(self)
        self._projection_ = _1Form_Projection(self)
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
                f"3dCSCG 1form FUNC do not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1form FUNC do not accept func {func_body.__class__}")

    @property
    def special(self):
        return self._special_

    @property
    def projection(self):
        """A wrapper of all projection methods."""
        return self._projection_

    def discretize(self, update_cochain=True, target = 'func', **kwargs):
        """
        Discretize the current function (a vector field:
        :class:`_3dCSCG.form.continuous.vector._3dCSCG_VectorField`) to cochain.
        It is actually a wrapper of multiple methods that discretize functions of different types (a vector
        field can be defined and represented in different ways in `python`, right?).
        :param target:
        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param kwargs:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if target == 'func':
            if self.TW.func.body.__class__.__name__ == '_3dCSCG_VectorField':

                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain, **kwargs)
                else:
                    raise NotImplementedError(f"3dCSCG 1-form cannot (target func) discretize _3dCSCG_VectorField of ftype={self.func.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 1-form can not discretize {self.TW.func.body.__class__}.')
        else:
            raise NotImplementedError(f"3dCSCG 1-form cannot discretize while targeting at {target}.")



    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain=True, quad_degree=None):
        if self.___DISCRETIZE_STANDARD_CACHE___ is None or \
            quad_degree != self.___DISCRETIZE_STANDARD_CACHE___['quadDegree']:
            self.___DISCRETIZE_STANDARD_CACHE___ = dict()

            xi, eta, sigma, edge_size_d_xi, quad_weights = \
                self.___PRIVATE_discretize_preparation___(d_='x', quad_degree=quad_degree)
            self.___DISCRETIZE_STANDARD_CACHE___['X'] = (xi, eta, sigma)

            xi, eta, sigma, edge_size_d_eta, quad_weights = \
                self.___PRIVATE_discretize_preparation___(d_='y', quad_degree=quad_degree)
            self.___DISCRETIZE_STANDARD_CACHE___['Y'] = (xi, eta, sigma)

            xi, eta, sigma, edge_size_d_sigma, quad_weights = \
                self.___PRIVATE_discretize_preparation___(d_='z', quad_degree=quad_degree)
            self.___DISCRETIZE_STANDARD_CACHE___['Z'] = (xi, eta, sigma)

            edge_size = (edge_size_d_xi, edge_size_d_eta, edge_size_d_sigma)
            self.___DISCRETIZE_STANDARD_CACHE___['edge'] = edge_size
            self.___DISCRETIZE_STANDARD_CACHE___['quad_weights'] = quad_weights
            self.___DISCRETIZE_STANDARD_CACHE___['quadDegree'] = quad_degree
        else:
            pass

        xi_x, eta_x, sigma_x = self.___DISCRETIZE_STANDARD_CACHE___['X']
        xi_y, eta_y, sigma_y = self.___DISCRETIZE_STANDARD_CACHE___['Y']
        xi_z, eta_z, sigma_z = self.___DISCRETIZE_STANDARD_CACHE___['Z']
        quad_weights = self.___DISCRETIZE_STANDARD_CACHE___['quad_weights']
        edge_size = self.___DISCRETIZE_STANDARD_CACHE___['edge']

        local_dx = dict()
        local_dy = dict()
        local_dz = dict()

        JXC, JYC, JZC = dict(), dict(), dict()
        for i in self.mesh.elements.indices:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark

            smctm = element.coordinate_transformation.mapping(xi_x, eta_x, sigma_x)
            if typeWr2Metric in JXC:
                J = JXC[typeWr2Metric]
            else:
                J = element.coordinate_transformation.Jacobian_matrix(xi_x, eta_x, sigma_x)
                if isinstance(typeWr2Metric, str):
                    JXC[typeWr2Metric] = J
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                u = self.func.body[0](*smctm)
                local_dx[i] = np.einsum(
                    'jk, j, k -> k', J[0][0]*u, quad_weights[0],
                    edge_size[0] * 0.5, optimize='greedy'
                )
            else:
                J = (J[0][0], J[1][0], J[2][0])
                u = self.func.body[0](*smctm)
                v = self.func.body[1](*smctm)
                w = self.func.body[2](*smctm)
                local_dx[i] = np.einsum(
                    'jk, j, k -> k', J[0]*u + J[1]*v + J[2]*w, quad_weights[0],
                    edge_size[0] * 0.5, optimize='greedy'
                )

            smctm = element.coordinate_transformation.mapping(xi_y, eta_y, sigma_y)
            if typeWr2Metric in JYC:
                J = JYC[typeWr2Metric]
            else:
                J = element.coordinate_transformation.Jacobian_matrix(xi_y, eta_y, sigma_y)
                if isinstance(typeWr2Metric, str):
                    JYC[typeWr2Metric] = J
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                v = self.func.body[1](*smctm)
                local_dy[i] = np.einsum(
                    'jk, j, k -> k', J[1][1]*v, quad_weights[1],
                    edge_size[1]*0.5, optimize='greedy'
                )
            else:
                J = (J[0][1], J[1][1], J[2][1])
                u = self.func.body[0](*smctm)
                v = self.func.body[1](*smctm)
                w = self.func.body[2](*smctm)
                local_dy[i] = np.einsum(
                    'jk, j, k -> k', J[0]*u + J[1]*v + J[2]*w, quad_weights[1],
                    edge_size[1]*0.5, optimize='greedy'
                )

            smctm = element.coordinate_transformation.mapping(xi_z, eta_z, sigma_z)
            if typeWr2Metric in JZC:
                J = JZC[typeWr2Metric]
            else:
                J = element.coordinate_transformation.Jacobian_matrix(xi_z, eta_z, sigma_z)
                if isinstance(typeWr2Metric, str):
                    JZC[typeWr2Metric] = J
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                w = self.func.body[2](*smctm)
                local_dz[i] = np.einsum(
                    'jk, j, k -> k', J[2][2]*w, quad_weights[2],
                    edge_size[2]*0.5, optimize='greedy'
                )
            else:
                J = (J[0][2], J[1][2], J[2][2])
                u = self.func.body[0](*smctm)
                v = self.func.body[1](*smctm)
                w = self.func.body[2](*smctm)
                local_dz[i] = np.einsum(
                    'jk, j, k -> k', J[0]*u + J[1]*v + J[2]*w, quad_weights[2],
                    edge_size[2]*0.5, optimize='greedy'
                )
        del JXC, JYC, JZC
        # isisKronecker? ...
        if not self.space.IS_Kronecker: raise NotImplementedError()
        # give it to cochain.local ...
        cochainLocal = dict()
        for i in self.mesh.elements.indices:
            cochainLocal[i] = np.hstack((local_dx[i], local_dy[i], local_dz[i]))
        if update_cochain:
            self.cochain.local = cochainLocal
        # ...
        return cochainLocal

    def ___PRIVATE_discretize_preparation___(self, d_='', quad_degree=None):
        p = [self.dqp[i] + 1 for i in range(self.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
        quad_num_nodes = [len(quad_nodes_i) for quad_nodes_i in quad_nodes]
        sbn0 = self.space.nodes[0]
        sbn1 = self.space.nodes[1]
        sbn2 = self.space.nodes[2]
        if d_ == 'x':
            a = sbn0[1:] - sbn0[:-1]
            a = a.ravel('F')
            b = (self.p[1] + 1) * (self.p[2] + 1)
            edge_size_x = np.tile(a, b)
            snbx = b * self.p[0]
            D = quad_nodes[0][:, np.newaxis].repeat(snbx, axis=1) + 1
            assert np.shape(D)[1] == len(edge_size_x)
            xi1 = D * edge_size_x / 2
            xi2 = np.tile(sbn0[:-1], b)
            xi = xi1 + xi2
            eta = np.tile(np.tile(sbn1[:, np.newaxis].repeat(quad_num_nodes[0], axis=1).T,
                                  (self.p[0], 1)).reshape((quad_num_nodes[0], self.p[0] * (self.p[1] + 1)),
                                                          order='F'),
                          (1, self.p[2] + 1))
            sigma = sbn2.repeat(self.p[0] * (self.p[1] + 1))[np.newaxis, :].repeat(
                quad_num_nodes[0], axis=0)
            return xi, eta, sigma, edge_size_x, quad_weights
        elif d_ == 'y':
            edge_size_y = np.tile(np.repeat((sbn1[1:] - sbn1[:-1]),
                                            self.p[0] + 1), self.p[2] + 1)
            xi = np.tile(sbn0, self.p[1] * (self.p[2] + 1))[np.newaxis, :].repeat(
                quad_num_nodes[1], axis=0)
            sn_by = self.NUM_basis_components[1]
            eta1 = (quad_nodes[1][:, np.newaxis].repeat(sn_by, axis=1) + 1) * edge_size_y / 2
            eta2 = np.tile(np.repeat(sbn1[:-1], (self.p[0] + 1)), (self.p[2] + 1))
            eta = eta1 + eta2
            sigma = sbn2.repeat(self.p[1] * (self.p[0] + 1))[np.newaxis, :].repeat(
                quad_num_nodes[1], axis=0)
            return xi, eta, sigma, edge_size_y, quad_weights
        elif d_ == 'z':
            edge_size_z = np.repeat((sbn2[1:] - sbn2[:-1]),
                                    self.p[0] + 1).repeat(self.p[1] + 1)
            xi = np.tile(sbn0, (self.p[1] + 1) * (self.p[2]))[np.newaxis, :].repeat(
                quad_num_nodes[2], axis=0)
            eta = np.tile(np.repeat(sbn1, (self.p[0] + 1)), self.p[2])[np.newaxis, :].repeat(
                quad_num_nodes[2], axis=0)
            sn_bz = self.NUM_basis_components[2]
            sigma1 = (quad_nodes[2][:, np.newaxis].repeat(sn_bz, axis=1) + 1) * edge_size_z / 2
            sigma2 = sbn2[:-1].repeat((self.p[0] + 1) * (self.p[1] + 1))
            sigma = sigma1 + sigma2
            return xi, eta, sigma, edge_size_z, quad_weights
        else:
            raise Exception()

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

        iJ = self.mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
        for i in INDICES:
            element = self.mesh.elements[i]
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            u = np.einsum('ij, i -> j', basis[0], self.cochain.local_('x')[i], optimize='optimal')
            v = np.einsum('ij, i -> j', basis[1], self.cochain.local_('y')[i], optimize='optimal')
            w = np.einsum('ij, i -> j', basis[2], self.cochain.local_('z')[i], optimize='optimal')
            value[i] = [None, None, None]
            typeWr2Metric = element.type_wrt_metric.mark
            iJi = iJ[i]
            if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
                value[i][0] = u*iJi[0][0]
                value[i][1] = v*iJi[1][1]
                value[i][2] = w*iJi[2][2]
            else:
                for j in range(3):
                    value[i][j] = u*iJi[0][j] + v*iJi[1][j] + w*iJi[2][j]
            if ravel:
                pass
            else:
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                # noinspection PyUnresolvedReferences
                value[i] = [value[i][j].reshape(shape, order='F') for j in range(3)]
        return xyz, value




    def ___OPERATORS_inner___(self, other, i, xietasigma, quad_weights, bfSelf, bfOther):
        """Note that here we only return a local matrix."""
        element = self.mesh.elements[i]
        mark = element.type_wrt_metric.mark
        J = element.coordinate_transformation.Jacobian_matrix(*xietasigma)
        sqrtg = element.coordinate_transformation.Jacobian(*xietasigma, J=J)
        iJ = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=J)
        g = element.coordinate_transformation.inverse_metric_matrix(*xietasigma, iJ=iJ)
        del J, iJ
        M00 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[0][0], bfOther[0], bfSelf[0])
        M11 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[1][1], bfOther[1], bfSelf[1])
        M22 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[2][2], bfOther[2], bfSelf[2])
        if isinstance(mark, str) and mark[:4] == 'Orth':
            M01 = None
            M02 = None
            M12 = None
            M10 = None
            M20 = None
            M21 = None
        else:
            M01 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[0][1], bfOther[0], bfSelf[1])
            M02 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[0][2], bfOther[0], bfSelf[2])
            M12 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[1][2], bfOther[1], bfSelf[2])
            if other is self:
                M10 = M01.T
                M20 = M02.T
                M21 = M12.T
            else:
                M10 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[1][0], bfOther[1], bfSelf[0])
                M20 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[2][0], bfOther[2], bfSelf[0])
                M21 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg*g[2][1], bfOther[2], bfSelf[1])
        Mi = spspa.bmat([(M00, M01, M02),
                         (M10, M11, M12),
                         (M20, M21, M22)], format='csc')
        return Mi

    @staticmethod
    def ___PRIVATE_inner_H1___(quad_weights, sqrt_g_g, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights*sqrt_g_g, bfO, bfS, optimize='optimal')
        return spspa.csc_matrix(M)


    def ___PRIVATE_basis_boundary_integral_over_region_side___(self,
        vector, region_name, side_name, quad_degree=None):
        """
        TODO: to be completed.

        :param vector:
        :param region_name:
        :param side_name:
        :param quad_degree:
        :return:
        """
        p = [self.dqp[i] + 3 for i in range(self.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad

        if side_name == 'F':
            quad_weights_2d = np.kron(quad_weights[1], quad_weights[0])


        v0, v1 = vector
        if isinstance(v0, (int, float)):
            pass
        else:
            raise NotImplementedError()
        if isinstance(v1, (int, float)):
            pass
        else:
            raise NotImplementedError()


        bfn_xi = self.space.basises[0].node_basis(x=quad_nodes[0])
        bfn_et = self.space.basises[1].node_basis(x=quad_nodes[1])


        bfe_xi = self.space.basises[0].edge_basis(x=quad_nodes[0])
        bfe_et = self.space.basises[1].edge_basis(x=quad_nodes[1])

        if side_name == 'F':
            B = (np.kron(bfn_et, bfe_xi), np.kron(bfe_et, bfn_xi))

        B0, B1 = B

        elements = self.mesh.DO.FIND_element_attach_to_region_side(region_name, side_name)
        elements = elements.ravel('F')
        ELEMENTS = list()
        for i in elements:
            if i in self.mesh.elements:
                ELEMENTS.append(i)


        GV = spspa.lil_matrix((1, self.GLOBAL_num_dofs))

        for i in ELEMENTS:
            dofs = self.numbering.___PRIVATE_DO_find_dofs_on_element_side___(i, side_name)
            mark = self.mesh.elements[i].type_wrt_metric.mark
            spacing = self.mesh.elements[i].spacing

            if mark[:4] == 'Orth':
                Lx, Ly, Lz = (spacing[0][1] - spacing[0][0]) * 2, \
                             (spacing[1][1] - spacing[1][0]) * 2, \
                             (spacing[2][1] - spacing[2][0]) * 2,

                # print(spacing, flush=True)
                Jx = Lx / 2
                Jy = Ly / 2
                Jz = Lz / 2

                if side_name in 'NS':
                    J0, J1 = Jy, Jz
                elif side_name in 'WE':
                    J0, J1 = Jx, Jz
                elif side_name in 'BF':
                    J0, J1 = Jx, Jy
                else:
                    raise Exception()

                if isinstance(v0, (int, float)):
                    # noinspection PyUnboundLocalVariable
                    cewti1 = np.einsum('ki, i -> k', B0*v1*J1, quad_weights_2d, optimize='greedy')
                if isinstance(v1, (int, float)):
                    cewti2 = - np.einsum('ki, i -> k', B1*v0*J0, quad_weights_2d, optimize='greedy')

                # noinspection PyUnboundLocalVariable
                GV[0, dofs] += np.concatenate((cewti1, cewti2))

            else:
                raise NotImplementedError(f"can not handle non-orthogonal elements yet!")

        GV = GlobalVector(GV.tocsr().T)
        return GV

    def ___OPERATORS_wedge___(self, f2, quad_degree=None):
        """"""
        assert self.ndim == f2.ndim and self.k + f2.k == 3, " <___STORAGE_OPERATORS_WEDGE___> "
        assert self.mesh == f2.mesh, "Meshes do not match."

        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], f2.dqp[i]])) for i in
                           range(3)]
        quad_nodes, _, quad_weights = self.space.DO_evaluate_quadrature(
            quad_degree)
        _, bfSelf = self.DO.evaluate_basis_at_meshgrid(
            *quad_nodes)
        _, bfOther = f2.DO.evaluate_basis_at_meshgrid(
            *quad_nodes)

        W00 = self.___PRIVATE_wedge_H1___(quad_weights, bfOther[0], bfSelf[0])
        W11 = self.___PRIVATE_wedge_H1___(quad_weights, bfOther[1], bfSelf[1])
        W22 = self.___PRIVATE_wedge_H1___(quad_weights, bfOther[2], bfSelf[2])
        W = spspa.bmat([(W00, None, None),
                        (None, W11, None),
                        (None, None, W22)], format='csc')
        return W


    @staticmethod
    def ___PRIVATE_wedge_H1___(quad_weights, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights, bfO, bfS, optimize='optimal')
        return spspa.csc_matrix(M)



class _1Form_Special(FrozenOnly):
    def __init__(self, _1sf):
        self._sf_ = _1sf
        self._vortex_detection_ = None
        self._freeze_self_()

    def cross_product(self, u, e, quad_degree=None):
        """
        We do ``(self X other, e)`` where ``self`` and ``other`` both are n-form, n be either 1 or 2.

        :return:
        """
        SCP_generator = ___3dCSCG_1Form_CrossProduct___(self._sf_, u, e, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, SCP_generator, 'no_cache')


    def ___PRIVATE_projected_into_2form_exactly___(self):
        """We project this 1form into a 2form exactly. Since it is a
        exact projection, we will use a space one degree higher than
        than the space of this 1form. The mesh will be the same mesh.
        """
        space = self._sf_.space
        mesh = self._sf_.mesh

        sp = space.p
        op = list()
        for i in sp: op.append(i+1)

        SPACE = space.__class__(op, None)

        f2 = _2Form(mesh, SPACE, is_hybrid=self._sf_.IS_hybrid,
            orientation=self._sf_.orientation,
            numbering_parameters=self._sf_.numbering._numbering_parameters_,
            name='Projected_2form_of_'+self._sf_.standard_properties.name)


        W21 = self._sf_.operators.wedge(f2)
        invM2 = f2.matrices.mass.inv

        lc1 = self._sf_.cochain.local

        lc2 = dict()
        for i in lc1:
            lc2[i] = invM2[i] @ W21[i] @ lc1[i]

        f2.cochain.local = lc2

        return f2

    @property
    def vortex_detection(self):
        if self._vortex_detection_ is None:
            self._vortex_detection_ = ___3dCSCG_1Form_Vortex_Detection___(self._sf_)
        return self._vortex_detection_

class ___3dCSCG_1Form_CrossProduct___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(\\omega \\times u, e)`.

    :param w1:
    :param u1:
    :param e1:
    :param quad_degree:
    """
    def __init__(self, w1, u1, e1, quad_degree=None):
        assert u1.ndim == w1.ndim == e1.ndim, " <_3dCSCG_1Form_CrossProduct> "
        assert u1.k == w1.k == e1.k == 1, " <_3dCSCG_1Form_CrossProduct> "
        assert u1.mesh == w1.mesh, "Meshes do not match."
        assert u1.mesh == e1.mesh, "Meshes do not match."
        self._mesh_ = u1.mesh
        self._u_ = u1
        self._w_ = w1
        self._e_ = e1

        if quad_degree is None:
            quad_degree = [int(np.max([u1.dqp[i], w1.dqp[i], e1.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = u1.space.DO_evaluate_quadrature(quad_degree)
        xietasigma, wbf = w1.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=True)
        _, ubf = u1.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        if e1 is u1:
            ebf = ubf
        else:
            _, ebf = e1.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        self._xietasigma_ = xietasigma
        self._qw_ = quad_weights
        self._wbf_ = wbf
        self._ubf_ = ubf
        self._ebf_ = ebf
        self._JM_ = self._mesh_.elements.coordinate_transformation.Jacobian_matrix(*xietasigma)
        self._iJ_ = self._mesh_.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=self._JM_)
        self.DO_reset_cache()
        self._freeze_self_()

    def DO_reset_cache(self):
        self._J_cache_ = dict()

    def _J_(self, i):
        element = self._mesh_.elements[i]
        typeWr2Metric = element.type_wrt_metric.mark
        if typeWr2Metric in self._J_cache_:
            return self._J_cache_[typeWr2Metric]
        else:
            JM = self._JM_[i]
            if isinstance(typeWr2Metric, str) and  typeWr2Metric[:4] == 'Orth':
                J00 = JM[1][1] * JM[2][2]
                J01 = 0
                J02 = 0
                J10 = 0
                J11 = JM[2][2] * JM[0][0]
                J12 = 0
                J20 = 0
                J21 = 0
                J22 = JM[0][0] * JM[1][1]
            else:
                J00 = JM[1][1] * JM[2][2] - JM[1][2] * JM[2][1]
                J01 = JM[2][1] * JM[0][2] - JM[2][2] * JM[0][1]
                J02 = JM[0][1] * JM[1][2] - JM[0][2] * JM[1][1]
                J10 = JM[1][2] * JM[2][0] - JM[1][0] * JM[2][2]
                J11 = JM[2][2] * JM[0][0] - JM[2][0] * JM[0][2]
                J12 = JM[0][2] * JM[1][0] - JM[0][0] * JM[1][2]
                J20 = JM[1][0] * JM[2][1] - JM[1][1] * JM[2][0]
                J21 = JM[2][0] * JM[0][1] - JM[2][1] * JM[0][0]
                J22 = JM[0][0] * JM[1][1] - JM[0][1] * JM[1][0]
            J = (J00, J01, J02, J10, J11, J12, J20, J21, J22)
            # cache it even for unique mesh cells (because we may use them multiple times when do temporal iterations.)
            self._J_cache_[typeWr2Metric] = J
            return J

    def __call__(self, i):
        typeWr2Metric = self._mesh_.elements[i].type_wrt_metric.mark

        u0, u1, u2 =  self._ubf_ # a
        w0, w1, w2 =  self._wbf_ # b; given
        e0, e1, e2 =  self._ebf_ # epsilon

        b0p = np.einsum('ij, i -> j', w0, self._w_.cochain.local_('x')[i], optimize='greedy')
        b1p = np.einsum('ij, i -> j', w1, self._w_.cochain.local_('y')[i], optimize='greedy')
        b2p = np.einsum('ij, i -> j', w2, self._w_.cochain.local_('z')[i], optimize='greedy')

        iJ = self._iJ_[i]
        J00, J01, J02, J10, J11, J12, J20, J21, J22 = self._J_(i)

        if isinstance(typeWr2Metric, str) and  typeWr2Metric[:4] == 'Orth':
            b0 = b0p * iJ[0][0]
            b1 = b1p * iJ[1][1]
            b2 = b2p * iJ[2][2]

            B01, B02 = -iJ[0][0]*b2, iJ[0][0]*b1
            B10, B12 = iJ[1][1]*b2, -iJ[1][1]*b0
            B20, B21 = -iJ[2][2]*b1, iJ[2][2]*b0

            m01 = B10*J00
            m02 = B20*J00
            m10 = B01*J11
            m12 = B21*J11
            m20 = B02*J22
            m21 = B12*J22

            # put `-` because a x b = - b x a
            M01 = - np.einsum('iw, jw, w -> ij', e0, u1, m01*self._qw_, optimize='greedy')
            M02 = - np.einsum('iw, jw, w -> ij', e0, u2, m02*self._qw_, optimize='greedy')

            M10 = - np.einsum('iw, jw, w -> ij', e1, u0, m10*self._qw_, optimize='greedy')
            M12 = - np.einsum('iw, jw, w -> ij', e1, u2, m12*self._qw_, optimize='greedy')

            M20 = - np.einsum('iw, jw, w -> ij', e2, u0, m20*self._qw_, optimize='greedy')
            M21 = - np.einsum('iw, jw, w -> ij', e2, u1, m21*self._qw_, optimize='greedy')

            M = ([None, spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
                 [spspa.csc_matrix(M10), None, spspa.csc_matrix(M12)],
                 [spspa.csc_matrix(M20), spspa.csc_matrix(M21), None])
        else:
            b0 = b0p * iJ[0][0] + b1p * iJ[1][0] + b2p * iJ[2][0]
            b1 = b0p * iJ[0][1] + b1p * iJ[1][1] + b2p * iJ[2][1]
            b2 = b0p * iJ[0][2] + b1p * iJ[1][2] + b2p * iJ[2][2]

            B00, B01, B02 = (iJ[0][1]*b2-iJ[0][2]*b1), (iJ[0][2]*b0-iJ[0][0]*b2), (iJ[0][0]*b1-iJ[0][1]*b0)
            B10, B11, B12 = (iJ[1][1]*b2-iJ[1][2]*b1), (iJ[1][2]*b0-iJ[1][0]*b2), (iJ[1][0]*b1-iJ[1][1]*b0)
            B20, B21, B22 = (iJ[2][1]*b2-iJ[2][2]*b1), (iJ[2][2]*b0-iJ[2][0]*b2), (iJ[2][0]*b1-iJ[2][1]*b0)

            m00 = B00*J00 + B01*J01 + B02*J02
            m01 = B10*J00 + B11*J01 + B12*J02
            m02 = B20*J00 + B21*J01 + B22*J02
            m10 = B00*J10 + B01*J11 + B02*J12
            m11 = B10*J10 + B11*J11 + B12*J12
            m12 = B20*J10 + B21*J11 + B22*J12
            m20 = B00*J20 + B01*J21 + B02*J22
            m21 = B10*J20 + B11*J21 + B12*J22
            m22 = B20*J20 + B21*J21 + B22*J22

            # put `-` because a x b = - b x a
            M00 = - np.einsum('iw, jw, w -> ij', e0, u0, m00*self._qw_, optimize='greedy')
            M01 = - np.einsum('iw, jw, w -> ij', e0, u1, m01*self._qw_, optimize='greedy')
            M02 = - np.einsum('iw, jw, w -> ij', e0, u2, m02*self._qw_, optimize='greedy')

            M10 = - np.einsum('iw, jw, w -> ij', e1, u0, m10*self._qw_, optimize='greedy')
            M11 = - np.einsum('iw, jw, w -> ij', e1, u1, m11*self._qw_, optimize='greedy')
            M12 = - np.einsum('iw, jw, w -> ij', e1, u2, m12*self._qw_, optimize='greedy')

            M20 = - np.einsum('iw, jw, w -> ij', e2, u0, m20*self._qw_, optimize='greedy')
            M21 = - np.einsum('iw, jw, w -> ij', e2, u1, m21*self._qw_, optimize='greedy')
            M22 = - np.einsum('iw, jw, w -> ij', e2, u2, m22*self._qw_, optimize='greedy')

            M = ([spspa.csc_matrix(M00), spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
                 [spspa.csc_matrix(M10), spspa.csc_matrix(M11), spspa.csc_matrix(M12)],
                 [spspa.csc_matrix(M20), spspa.csc_matrix(M21), spspa.csc_matrix(M22)])
        MW = spspa.bmat(M, format='csc')

        return MW

class ___3dCSCG_1Form_Vortex_Detection___(FrozenOnly):
    """A wrapper of all vortex detection methods. So, we consider this 1 form as
    a variable of a flow field."""
    def __init__(self, _1sf):
        self._sf_ = _1sf
        self._freeze_self_()

    def ___PRIVATE_generate_gradient_tensor_at___(self, xi, eta, sigma):
        """We compute the gradient tensor of this 1form. To do so, we first project
        this 1-form into a vector of 3 standard 0-forms which represent the three
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






class _1Form_Projection(FrozenOnly):
    """A wrapper of all projection methods."""
    def __init__(self,_1sf):
        self._sf_ = _1sf
        self._to_ = ___3dCSCG_1Form_Project_To___(self._sf_)
        self._van_ = ___3dCSCG_1Form_Project_Van___(self._sf_)
        self._freeze_self_()

    @property
    def to(self):
        return self._to_

    @property
    def van(self):
        return self._van_


class ___3dCSCG_1Form_Project_To___(FrozenOnly):
    """A wrapper of all projection into methods"""
    def __init__(self,_1sf):
        self._sf_ = _1sf
        self._freeze_self_()

    def standard_2form(self):
        """Project this 1form into a 2form exactly."""
        return self._sf_.special.___PRIVATE_projected_into_2form_exactly___()

    def vector_of_3_standard_0forms(self):
        """project this 1form into a tuple of three 0forms. Each 0form stands for
        a component of the 1form (as a vector).

        Since the 0forms with the same space will be of higher degree, so we do not
        need to make a new space of higher degree. And we of course use the same mesh.
        Thus, both the mesh and space will the same as the those of this 1form."""
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


class ___3dCSCG_1Form_Project_Van___(FrozenOnly):
    """A wrapper of all projection into methods"""
    def __init__(self,_1sf):
        self._sf_ = _1sf
        self._freeze_self_()






if __name__ == '__main__':
    # mpiexec -n 4 python _3dCSCG\form\standard\_1_form.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([8,8,8])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    velocity = FC('vector', (u,v,w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)


    f1 = FC('1-f', is_hybrid=False)


    f1.TW.func.DO.set_func_body_as(velocity)
    f1.TW.current_time = 0
    f1.TW.___DO_push_all_to_instant___()
    f1.discretize()

    VD = f1.special.vortex_detection
    SO = VD.Q_and_lambda2([-1, 0, 0.5, 1], [-1, 1], [-1, -0.5, 0, 0.5, 1])