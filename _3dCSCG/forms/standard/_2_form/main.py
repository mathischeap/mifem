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
from _3dCSCG.forms.standard._2_form.discretize.main import _3dCSCG_Discretize
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
        self._discretize_ = _3dCSCG_Discretize(self)
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
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

    @property
    def discretize(self):
        return self._discretize_






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
            u = np.einsum('ij, i -> j', basis[0], self.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
            v = np.einsum('ij, i -> j', basis[1], self.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')
            w = np.einsum('ij, i -> j', basis[2], self.cochain.___PRIVATE_local_on_axis___('z', i), optimize='greedy')
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
                    _0u = iJi[1][1] * iJi[2][2] - iJi[1][2] * iJi[2][1]
                    _1v = iJi[2][2] * iJi[0][0] - iJi[2][0] * iJi[0][2]
                    _2w = iJi[0][0] * iJi[1][1] - iJi[0][1] * iJi[1][0]
                    rm00 = np.einsum('ji, j -> ji', b0, _0u, optimize='greedy')
                    rm11 = np.einsum('ji, j -> ji', b1, _1v, optimize='greedy')
                    rm22 = np.einsum('ji, j -> ji', b2, _2w, optimize='greedy')

                    if typeWr2Metric[:4] == 'Orth':
                        RM_i_ = ( np.hstack((rm00, OO01, OO02)),
                                  np.hstack((OO10, rm11, OO12)),
                                  np.hstack((OO20, OO21, rm22)) )
                    else:
                        _0v = iJi[2][1] * iJi[0][2] - iJi[2][2] * iJi[0][1]
                        _0w = iJi[0][1] * iJi[1][2] - iJi[0][2] * iJi[1][1]
                        _1u = iJi[1][2] * iJi[2][0] - iJi[1][0] * iJi[2][2]
                        _1w = iJi[0][2] * iJi[1][0] - iJi[0][0] * iJi[1][2]
                        _2u = iJi[1][0] * iJi[2][1] - iJi[1][1] * iJi[2][0]
                        _2v = iJi[2][0] * iJi[0][1] - iJi[2][1] * iJi[0][0]
                        rm01 = np.einsum('ji, j -> ji', b1, _0v, optimize='greedy')
                        rm02 = np.einsum('ji, j -> ji', b2, _0w, optimize='greedy')
                        rm10 = np.einsum('ji, j -> ji', b0, _1u, optimize='greedy')
                        rm12 = np.einsum('ji, j -> ji', b2, _1w, optimize='greedy')
                        rm20 = np.einsum('ji, j -> ji', b0, _2u, optimize='greedy')
                        rm21 = np.einsum('ji, j -> ji', b1, _2v, optimize='greedy')
                        RM_i_ = ( np.hstack((rm00, rm01, rm02)),
                                  np.hstack((rm10, rm11, rm12)),
                                  np.hstack((rm20, rm21, rm22)))

                    CACHE[typeWr2Metric] = RM_i_
                    RM[i] = RM_i_
            else:
                iJi = iJ[i]
                _0u = iJi[1][1] * iJi[2][2] - iJi[1][2] * iJi[2][1]
                _0v = iJi[2][1] * iJi[0][2] - iJi[2][2] * iJi[0][1]
                _0w = iJi[0][1] * iJi[1][2] - iJi[0][2] * iJi[1][1]
                _1u = iJi[1][2] * iJi[2][0] - iJi[1][0] * iJi[2][2]
                _1v = iJi[2][2] * iJi[0][0] - iJi[2][0] * iJi[0][2]
                _1w = iJi[0][2] * iJi[1][0] - iJi[0][0] * iJi[1][2]
                _2u = iJi[1][0] * iJi[2][1] - iJi[1][1] * iJi[2][0]
                _2v = iJi[2][0] * iJi[0][1] - iJi[2][1] * iJi[0][0]
                _2w = iJi[0][0] * iJi[1][1] - iJi[0][1] * iJi[1][0]


                rm00 = np.einsum('ji, j -> ji', b0, _0u, optimize='greedy')
                rm01 = np.einsum('ji, j -> ji', b1, _0v, optimize='greedy')
                rm02 = np.einsum('ji, j -> ji', b2, _0w, optimize='greedy')

                rm10 = np.einsum('ji, j -> ji', b0, _1u, optimize='greedy')
                rm11 = np.einsum('ji, j -> ji', b1, _1v, optimize='greedy')
                rm12 = np.einsum('ji, j -> ji', b2, _1w, optimize='greedy')

                rm20 = np.einsum('ji, j -> ji', b0, _2u, optimize='greedy')
                rm21 = np.einsum('ji, j -> ji', b1, _2v, optimize='greedy')
                rm22 = np.einsum('ji, j -> ji', b2, _2w, optimize='greedy')

                RM[i] = ( np.hstack((rm00, rm01, rm02)),
                          np.hstack((rm10, rm11, rm12)),
                          np.hstack((rm20, rm21, rm22)))

        return RM


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
    def ___PRIVATE_inner_H1___(quad_weights, sqrt_g, g, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights * sqrt_g * g, bfO, bfS, optimize='optimal')
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
