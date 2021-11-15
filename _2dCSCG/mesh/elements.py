# -*- coding: utf-8 -*-

from SCREWS.frozen import FrozenOnly
from SCREWS.decorators import memoize1, accepts
from SCREWS.quadrature import Quadrature
from root.config import *



class _2dCSCG_Mesh_Elements(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = dict()
        self._ct_ = _2dCSCG_Mesh_Elements_CT(self)
        for i in self.indices:
            self._elements_[i] = _2dCSCG_Mesh_Element(self, i)
        self.___PRIVATE_parse_elements_type_wrt_metric___()
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        """"""
        self.coordinate_transformation.RESET_cache()

    def ___PRIVATE_parse_elements_type_wrt_metric___(self):
        counter: dict = dict()
        self._multi_elements_metric_: dict = dict()
        for i in self:
            ei = self[i]
            mki = ei.type_wrt_metric.mark
            if mki in counter:
                ei._type_wrt_metric_ = self[counter[mki]]._type_wrt_metric_
                if mki in self._multi_elements_metric_:
                    self._multi_elements_metric_[mki] += 1
                else:
                    self._multi_elements_metric_[mki] = 2
            else:
                counter[mki] = i

    @property
    def GLOBAL_num(self):
        return self._mesh_._num_total_elements_

    @property
    def in_regions(self):
        """The regions elements in this core are in."""
        return self._mesh_._elements_in_regions_

    @property
    def indices(self):
        return self._mesh_._element_indices_

    @property
    def num(self):
        return self._mesh_._num_local_elements_

    @property
    def map(self):
        return self._mesh_.___element_map___

    @property
    def layout(self):
        return self._mesh_._element_layout_

    @property
    def spacing(self):
        return self._mesh_._element_spacing_

    @property
    def ratio(self):
        return self._mesh_._element_ratio_

    @property
    def coordinate_transformation(self):
        return self._ct_


    def __getitem__(self, i):
        return self._elements_[i]

    def __contains__(self, i):
        return i in self.indices

    def __iter__(self):
        for i in self.indices:
            yield i

    def __len__(self):
        return len(self.indices)

    def __eq__(self, other):
        assert other.__class__.__name__ == '_2dCSCG_Mesh_Elements', \
            'I can only equal to a _3dCSCG_Mesh_Elements object.'
        # when mesh is equal, elements must be
        mesh_judge = self._mesh_ == other._mesh_
        indices_judge = self.indices == other.indices
        judge = mesh_judge and indices_judge
        return cOmm.allreduce(judge, op=MPI.LAND)


    def ___PRIVATE_elementwise_cache_metric_key___(self, i):
        mark = self[i].type_wrt_metric.mark
        if not isinstance(mark, str): mark = str(mark)
        return mark

    @accepts('self', int)
    def ___DO_find_slave_of_element___(self, i: int) -> int:
        return self._mesh_.DO.FIND_slave_of_element(i)





class _2dCSCG_Mesh_Element(FrozenOnly):
    """"""
    def __init__(self, elements, i):
        self._elements_ = elements
        self._mesh_ = elements._mesh_
        self._i_ = i
        self._type_wrt_metric_ = None
        self._in_region_ = self._mesh_.DO.FIND_region_name_of_element(self.i)
        self._ct_ = _2dCSCG_Mesh_ECT(self)
        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def position(self):
        return self._elements_.map[self.i]

    @property
    def in_region(self):
        return self._in_region_

    @property
    def spacing(self):
        region, localRegionIndices = self._mesh_.DO.FIND_region_name_and_local_indices_of_element(self.i)
        elementsSpacing = self._elements_.spacing[region]
        _spacing_ = np.zeros((2,2))
        for i in range(2):
            _spacing_[i, 0] = elementsSpacing[i][localRegionIndices[i]]
            _spacing_[i, 1] = elementsSpacing[i][localRegionIndices[i]+1]
        return _spacing_

    @property
    def type_wrt_metric(self):
        if self._type_wrt_metric_ is None:
            region, _ = self._mesh_.DO.FIND_region_name_and_local_indices_of_element(self.i)
            self._type_wrt_metric_ = \
                self._mesh_.domain.regions[region].type_wrt_metric.___CLASSIFY_ELEMENT_of_spacing___(
                    self.spacing)
        return self._type_wrt_metric_

    @property
    def coordinate_transformation(self):
        return self._ct_





class _2dCSCG_Mesh_ECT(FrozenOnly):
    def __init__(self, element):
        self._element_ = element
        self._mesh_ = self._element_._elements_._mesh_
        self._region_ = self._mesh_.domain.regions[self._element_.in_region]
        self._origin_ = None
        self._delta_ = None
        self._freeze_self_()

    @property
    def origin(self):
        if self._origin_ is None:
            in_region, local_indices = self._mesh_.DO.FIND_region_name_and_local_indices_of_element(
                self._element_.i)
            self._origin_, self._delta_ = \
                self._mesh_.DO.FIND_reference_origin_and_size_of_element_of_given_local_indices(
                in_region, local_indices)
        return self._origin_

    @property
    def delta(self):
        if self._delta_ is None:
            in_region, local_indices = self._mesh_.DO.FIND_region_name_and_local_indices_of_element(
                self._element_.i)
            self._origin_, self._delta_ = \
                self._mesh_.DO.FIND_reference_origin_and_size_of_element_of_given_local_indices(
                in_region, local_indices)
        return self._delta_


    def mapping(self, *evaluationPoints):
        XY = self._region_.interpolation(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)])
        return XY

    def X(self, *evaluationPoints):
        X = self._region_.interpolation.mapping_X(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)])
        return X
    def Y(self, *evaluationPoints):
        Y = self._region_.interpolation.mapping_Y(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)])
        return Y


    def Jacobian_matrix(self, *evaluationPoints):
        mark = self._element_.type_wrt_metric.mark
        if isinstance(mark, str) and mark[:4] == 'Orth':
            xyz_xietasigma = [[0, 0], [0, 0]]
            rs = [(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)]
            J00 = self._region_.interpolation.Jacobian_Xr(*rs)
            J11 = self._region_.interpolation.Jacobian_Ys(*rs)
            xyz_xietasigma[0][0] = J00 * (self.delta[0] / 2)
            xyz_xietasigma[1][1] = J11 * (self.delta[1] / 2)
        else:
            xyz_xietasigma = [[np.zeros(np.shape(evaluationPoints[j])) for _ in range(2)] for j in range(2)]
            xyz_rst = self._region_.interpolation.Jacobian_matrix(
                *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)])
            for j in range(2):
                for l in range(2):
                    xyz_xietasigma[j][l] = xyz_rst[j][l] * (self.delta[l] / 2)
        return xyz_xietasigma



    def J00(self, *evaluationPoints):
        Xr = self._region_.interpolation.Jacobian_Xr(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)])
        return Xr * self.delta[0] / 2
    def J01(self, *evaluationPoints):
        Xs = self._region_.interpolation.Jacobian_Xs(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)])
        return Xs * self.delta[1] / 2

    def J10(self, *evaluationPoints):
        Yr = self._region_.interpolation.Jacobian_Yr(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)])
        return Yr * self.delta[0] / 2
    def J11(self, *evaluationPoints):
        Ys = self._region_.interpolation.Jacobian_Ys(
            *[(evaluationPoints[j] + 1) * 0.5 * self.delta[j] + self.origin[j] for j in range(2)])
        return Ys * self.delta[1] / 2

    def J0_(self, *evaluationPoints):
        return self.J00(*evaluationPoints), self.J01(*evaluationPoints)
    def J1_(self, *evaluationPoints):
        return self.J10(*evaluationPoints), self.J11(*evaluationPoints)






    def Jacobian(self, *evaluationPoints, J=None):
        """Determinant of the Jacobian matrix."""
        if J is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        return J[0][0]*J[1][1] - J[0][1]*J[1][0]

    def metric(self, *evaluationPoints, detJ=None):
        """
        The metric ``g:= det(G):=(det(J))**2``. Since our Jacobian and inverse of Jacobian are both square,
        we know that the metric ``g`` is equal to square of ``det(J)``. ``g = (det(J))**2`` is due to the
        fact that the Jacobian matrix is square. The definition of ``g`` usually is given
        as ``g:= det(G)`` where ``G`` is the metric matrix, or metric tensor.
        """
        if detJ is None:
            detJ = self.Jacobian(*evaluationPoints)
        return detJ ** 2



    def inverse_Jacobian_matrix(self, *evaluationPoints, J=None):
        """The inverse Jacobian matrix. """
        if J is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        Jacobian = J[0][0]*J[1][1] - J[0][1]*J[1][0]
        reciprocalJacobian = 1 / Jacobian
        del Jacobian
        iJ00 = + reciprocalJacobian * J[1][1]
        iJ01 = - reciprocalJacobian * J[0][1]
        iJ10 = - reciprocalJacobian * J[1][0]
        iJ11 = + reciprocalJacobian * J[0][0]
        return [[iJ00, iJ01],
                [iJ10, iJ11]]

    def inverse_Jacobian(self, *evaluationPoints, iJ=None):
        """Determinant of the inverse Jacobian matrix. """
        if iJ is None:
            iJ = self.inverse_Jacobian_matrix(*evaluationPoints)
        return iJ[0][0]*iJ[1][1] - iJ[0][1]*iJ[1][0]



    def metric_matrix(self, *evaluationPoints, J=None):
        """
        Also called metric tensor. Let J be the Jacobian matrix. The ``metricMatrix`` is
        denoted by G, G := J^T.dot(J). And the metric is ``g := (det(J))**2 or g := det(G).``
        Which means for a square Jacobian matrix, the metric turns out to be the square of the
        determinant of the Jacobian matrix.

        The entries of G is normally denoted as g_{i,j}.
        """
        if J is None:
            J = self.Jacobian_matrix(*evaluationPoints)
        G = [[None for _ in range(2)] for __ in range(2)]
        for i in range(2):
            for j in range(i, 2):
                # noinspection PyTypeChecker
                G[i][j] = J[0][i] * J[0][j]
                for l in range(1, 2):
                    G[i][j] += J[l][i] * J[l][j]
                if i != j:
                    G[j][i] = G[i][j]
        return G

    def inverse_metric_matrix(self, *evaluationPoints, iJ=None):
        """
        The ``inverseMetricMatrix`` is the metric matrix of the inverse Jacobian matrix
        or the metric of the inverse mapping. It is usually denoted as G^{-1}.

        The entries of G^{-1} is normally denoted as g^{i,j}.
        """
        if iJ is None:
            iJ = self.inverse_Jacobian_matrix(*evaluationPoints)
        iG = [[None for _ in range(2)] for __ in range(2)]
        for i in range(2):
            for j in range(i, 2):
                # noinspection PyTypeChecker
                iG[i][j] = iJ[i][0] * iJ[j][0]
                for l in range(1, 2):
                    iG[i][j] += iJ[i][l] * iJ[j][l]
                if i != j:
                    iG[j][i] = iG[i][j]
        return iG





class _2dCSCG_Mesh_Elements_CT(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self._ctq_1d_ = _2dCSCG_ECT_1d_QUAD(self)
        self._ctq_2d_ = _2dCSCG_ECT_2d_QUAD(self)

    @property
    def QUAD_1d(self):
        """Evaluating the coordinate transformation from quadrature. Results are in 1d array (do ravel)."""
        return self._ctq_1d_

    @property
    def QUAD_2d(self):
        """Evaluating the coordinate transformation from quadrature. Results are in 2d array."""
        return self._ctq_2d_


    def mapping(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta)
    def X(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta)
    def Y(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta)

    def Jacobian_matrix(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta)
    def J00(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta)
    def J01(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta)
    def J10(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta)
    def J11(self, xi, eta):
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta)

    def Jacobian(self, xi, eta, J=None):
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, intermediateData=J)
    def metric(self, xi, eta, detJ=None):
        """g := det(G) = Jacobian ** 2."""
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, intermediateData=detJ)
    def metric_matrix(self, xi, eta, J=None):
        """G, g_{i,j}."""
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, intermediateData=J)

    def inverse_Jacobian_matrix(self, xi, eta, J=None):
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, intermediateData=J)
    def inverse_Jacobian(self, xi, eta, iJ=None):
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, intermediateData=iJ)
    def inverse_metric_matrix(self, xi, eta, iJ=None):
        """G^-1, g^{i,j}."""
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, intermediateData=iJ)






class ElementsCTValuesCache(FrozenOnly):
    def __init__(self, elements, CTT, xi, eta, intermediateData=None):
        self._elements_ = elements
        self._CTT_ = CTT
        self._xi_eta_sigma_ = (xi, eta)
        if intermediateData is not None:
            assert intermediateData.__class__.__name__ == 'ElementsCTValuesCache', \
                "intermediateData can only be an ElementsCTValues object."
            inter_xi_eta_sigma = intermediateData._xi_eta_sigma_
            for i, com in enumerate(inter_xi_eta_sigma):
                assert self._xi_eta_sigma_[i] is com, \
                    "intermediateData xi, eta sigma must be the same object."
            if CTT == 'Jacobian':
                assert intermediateData._CTT_ == 'Jacobian_matrix'
            elif CTT == 'metric':
                assert intermediateData._CTT_ == 'Jacobian'
            elif CTT == 'metric_matrix':
                assert intermediateData._CTT_ == 'Jacobian_matrix'
            elif CTT == 'inverse_Jacobian_matrix':
                assert intermediateData._CTT_ == 'Jacobian_matrix'
            elif CTT == 'inverse_Jacobian':
                assert intermediateData._CTT_ == 'inverse_Jacobian_matrix'
            elif CTT == 'inverse_metric_matrix':
                assert intermediateData._CTT_ == 'inverse_Jacobian_matrix'
            else:
                raise Exception(f"{CTT} ask for no intermediateData or provided wrong intermediateData.")
        else:
            pass
        self._intermediateData_ = intermediateData
        self._multi_elements_metric_ = self._elements_._multi_elements_metric_
        self._cache_ = dict()
        self._freeze_self_()

    def __getitem__(self, i):
        element = self._elements_[i]
        type_wrt_metric = element.type_wrt_metric.mark
        if self._CTT_ in ('mapping', 'X', 'Y'):
            return getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
        elif self._CTT_ in ('Jacobian_matrix', 'J00', 'J01', 'J10', 'J11'):
            if type_wrt_metric in self._cache_:
                return self._cache_[type_wrt_metric]
            else:
                JM = getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
                if isinstance(type_wrt_metric, str) and \
                    type_wrt_metric in self._multi_elements_metric_ and \
                    type_wrt_metric not in self._cache_ and \
                    self._multi_elements_metric_[type_wrt_metric] >= caChe_factor:
                    self._cache_[type_wrt_metric] = JM
                return JM
        else:
            if type_wrt_metric in self._cache_:
                return self._cache_[type_wrt_metric]
            else:
                if self._intermediateData_ is None:
                    result = getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_)
                else:
                    result = getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_, self._intermediateData_[i])
                if isinstance(type_wrt_metric, str) and \
                    type_wrt_metric in self._multi_elements_metric_ and \
                    type_wrt_metric not in self._cache_ and \
                    self._multi_elements_metric_[type_wrt_metric] >= caChe_factor:
                    # here we have very strict cache rule.
                    self._cache_[type_wrt_metric] = result
                return result

    def __len__(self):
        return len(self._elements_)

    def __contains__(self, item):
        return item in self._elements_

    def __iter__(self):
        for i in self._elements_:
            yield i





class _2dCSCG_ECT_1d_QUAD(FrozenOnly):
    def __init__(self, ect):
        self._elements_ = ect._elements_
        self._freeze_self_()

    @staticmethod
    def ___compute_xietasigma___(quad_degree, quad_type):
        _Quadrature_ = Quadrature(quad_degree, category=quad_type)
        quad_nodes = _Quadrature_.quad[0]
        xi, eta = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        return xi, eta

    @memoize1
    def mapping(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta)

    @memoize1
    def X(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta)
    @memoize1
    def Y(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta)

    @memoize1
    def Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta)
    @memoize1
    def J00(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta)
    @memoize1
    def J01(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta)
    @memoize1
    def J10(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta)
    @memoize1
    def J11(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta)

    @memoize1
    def Jacobian(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, intermediateData=None)
    @memoize1
    def metric(self, quad_degree, quad_type):
        """g := det(G) = Jacobian ** 2."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, intermediateData=None)
    @memoize1
    def metric_matrix(self, quad_degree, quad_type):
        """G, g_{i,j}."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, intermediateData=None)

    @memoize1
    def inverse_Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, intermediateData=None)
    @memoize1
    def inverse_Jacobian(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, intermediateData=None)
    @memoize1
    def inverse_metric_matrix(self, quad_degree, quad_type):
        """G^-1, g^{i,j}."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, intermediateData=None)



class _2dCSCG_ECT_2d_QUAD(FrozenOnly):
    def __init__(self, ect):
        self._elements_ = ect._elements_
        self._freeze_self_()

    @staticmethod
    def ___compute_xietasigma___(quad_degree, quad_type):
        _Quadrature_ = Quadrature(quad_degree, category=quad_type)
        quad_nodes = _Quadrature_.quad[0]
        xi, eta = np.meshgrid(*quad_nodes, indexing='ij')
        return xi, eta

    @memoize1
    def mapping(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta)

    @memoize1
    def X(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta)
    @memoize1
    def Y(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta)

    @memoize1
    def Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta)
    @memoize1
    def J00(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta)
    @memoize1
    def J01(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta)
    @memoize1
    def J10(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta)
    @memoize1
    def J11(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta)

    @memoize1
    def Jacobian(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, intermediateData=None)
    @memoize1
    def metric(self, quad_degree, quad_type):
        """g := det(G) = Jacobian ** 2."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, intermediateData=None)
    @memoize1
    def metric_matrix(self, quad_degree, quad_type):
        """G, g_{i,j}."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, intermediateData=None)

    @memoize1
    def inverse_Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, intermediateData=None)
    @memoize1
    def inverse_Jacobian(self, quad_degree, quad_type):
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, intermediateData=None)
    @memoize1
    def inverse_metric_matrix(self, quad_degree, quad_type):
        """G^-1, g^{i,j}."""
        xi, eta = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, intermediateData=None)