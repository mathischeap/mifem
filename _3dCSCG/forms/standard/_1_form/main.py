# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('../')
from root.config.main import *
from scipy import sparse as spspa
from screws.quadrature import Quadrature
from _3dCSCG.forms.standard._1_form.discretize.main import _3dCSCG_Discretize
from _3dCSCG.forms.standard.base.main import _3dCSCG_Standard_Form
from tools.linear_algebra.data_structures.global_matrix.main import GlobalVector
from _3dCSCG.forms.standard._1_form.project.main import _1Form_Projection
from _3dCSCG.forms.standard._1_form.special.main import _1Form_Special

class _3dCSCG_1Form(_3dCSCG_Standard_Form):
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
        self.___PRIVATE_reset_cache___()
        self._discretize_ = _3dCSCG_Discretize(self)
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
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
                ri = self.mesh.do.FIND_region_name_of_element(i)
                if ri in regions:
                    INDICES.append(i)

        iJ = self.mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
        for i in INDICES:
            element = self.mesh.elements[i]
            xyz[i] = element.coordinate_transformation.mapping(*xietasigma)
            u = np.einsum('ij, i -> j', basis[0], self.cochain.___PRIVATE_local_on_axis___('x', i), optimize='optimal')
            v = np.einsum('ij, i -> j', basis[1], self.cochain.___PRIVATE_local_on_axis___('y', i), optimize='optimal')
            w = np.einsum('ij, i -> j', basis[2], self.cochain.___PRIVATE_local_on_axis___('z', i), optimize='optimal')
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
                # noinspection PyUnresolvedReferences
                xyz[i] = [xyz[i][j].reshape(shape, order='F') for j in range(3)]
                # noinspection PyUnresolvedReferences
                value[i] = [value[i][j].reshape(shape, order='F') for j in range(3)]
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
        iJ = self.mesh.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma)
        b0, b1, b2 = basis
        b0 = b0.T
        b1 = b1.T
        b2 = b2.T
        OO01 = 0 * b1
        OO02 = 0 * b2
        OO10 = 0 * b0
        OO12 = 0 * b2
        OO20 = 0 * b0
        OO21 = 0 * b1
        CACHE = dict()
        RM = dict()
        for i in self.mesh.elements:
            element = self.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in CACHE:
                    RM[i] = CACHE[typeWr2Metric]
                else:
                    iJi = iJ[i]
                    rm00 = np.einsum('ji, j -> ji', b0, iJi[0][0], optimize='greedy')
                    rm11 = np.einsum('ji, j -> ji', b1, iJi[1][1], optimize='greedy')
                    rm22 = np.einsum('ji, j -> ji', b2, iJi[2][2], optimize='greedy')

                    if typeWr2Metric[:4] == 'Orth':
                        RM_i_ = ( np.hstack((rm00, OO01, OO02)),
                                  np.hstack((OO10, rm11, OO12)),
                                  np.hstack((OO20, OO21, rm22)) )
                    else:
                        rm01 = np.einsum('ji, j -> ji', b1, iJi[1][0], optimize='greedy')
                        rm02 = np.einsum('ji, j -> ji', b2, iJi[2][0], optimize='greedy')
                        rm10 = np.einsum('ji, j -> ji', b0, iJi[0][1], optimize='greedy')
                        rm12 = np.einsum('ji, j -> ji', b2, iJi[2][1], optimize='greedy')
                        rm20 = np.einsum('ji, j -> ji', b0, iJi[0][2], optimize='greedy')
                        rm21 = np.einsum('ji, j -> ji', b1, iJi[1][2], optimize='greedy')
                        RM_i_ = ( np.hstack((rm00, rm01, rm02)),
                                  np.hstack((rm10, rm11, rm12)),
                                  np.hstack((rm20, rm21, rm22)))

                    CACHE[typeWr2Metric] = RM_i_
                    RM[i] = RM_i_
            else:
                iJi = iJ[i]

                rm00 = np.einsum('ji, j -> ji', b0, iJi[0][0], optimize='greedy')
                rm01 = np.einsum('ji, j -> ji', b1, iJi[1][0], optimize='greedy')
                rm02 = np.einsum('ji, j -> ji', b2, iJi[2][0], optimize='greedy')

                rm10 = np.einsum('ji, j -> ji', b0, iJi[0][1], optimize='greedy')
                rm11 = np.einsum('ji, j -> ji', b1, iJi[1][1], optimize='greedy')
                rm12 = np.einsum('ji, j -> ji', b2, iJi[2][1], optimize='greedy')

                rm20 = np.einsum('ji, j -> ji', b0, iJi[0][2], optimize='greedy')
                rm21 = np.einsum('ji, j -> ji', b1, iJi[1][2], optimize='greedy')
                rm22 = np.einsum('ji, j -> ji', b2, iJi[2][2], optimize='greedy')

                RM[i] = ( np.hstack((rm00, rm01, rm02)),
                          np.hstack((rm10, rm11, rm12)),
                          np.hstack((rm20, rm21, rm22)))

        return RM




    def ___PRIVATE_operator_inner___(self, other, i, xietasigma, quad_weights, bfSelf, bfOther):
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

        elements = self.mesh.do.FIND_element_attach_to_region_side(region_name, side_name)
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

    def ___PRIVATE_operator_wedge___(self, f2, quad_degree=None):
        """"""
        assert self.ndim == f2.ndim and self.k + f2.k == 3, " <___STORAGE_OPERATORS_WEDGE___> "
        assert self.mesh == f2.mesh, "Meshes do not match."

        if quad_degree is None:
            quad_degree = [int(np.max([self.dqp[i], f2.dqp[i]])) for i in
                           range(3)]
        quad_nodes, _, quad_weights = self.space.___PRIVATE_do_evaluate_quadrature___(
            quad_degree)
        _, bfSelf = self.do.evaluate_basis_at_meshgrid(
            *quad_nodes)
        _, bfOther = f2.do.evaluate_basis_at_meshgrid(
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


    f1.TW.func.do.set_func_body_as(velocity)
    f1.TW.current_time = 0
    f1.TW.___DO_push_all_to_instant___()
    f1.discretize()

    VD = f1.special.vortex_detection
    SO = VD.Q_and_lambda2([-1, 0, 0.5, 1], [-1, 1], [-1, -0.5, 0, 0.5, 1])