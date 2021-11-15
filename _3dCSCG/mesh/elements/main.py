# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')

from SCREWS.frozen import FrozenOnly
from SCREWS.decorators import accepts
from SCREWS.quadrature import Quadrature
from root.config import *
from _3dCSCG.mesh.elements.element import _3dCSCG_Mesh_Element

import matplotlib.pyplot as plt


class _3dCSCG_Mesh_Elements(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = dict()
        self._DO_ = None
        self._ct_ = _3dCSCG_Mesh_Elements_CT(self)
        for i in self.indices:
            self._elements_[i] = _3dCSCG_Mesh_Element(self, i)
        self.___PRIVATE_parse_elements_type_wrt_metric___()
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self.coordinate_transformation.RESET_cache()
        self._quality_ = None

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
        """The total number of elements in all cores."""
        return self._mesh_._num_total_elements_

    @property
    def in_regions(self):
        """The elements in this core are in these regions."""
        return self._mesh_._elements_in_regions_

    @property
    def indices(self):
        """The local elements are of these indices."""
        return self._mesh_._element_indices_

    @property
    def num(self):
        """How many elements in this local core."""
        return self._mesh_._num_local_elements_

    @property
    def map(self):
        """The element map for local elements."""
        return self._mesh_.___element_map___

    @property
    def layout(self):
        """The element layout. A dict."""
        return self._mesh_._element_layout_

    @property
    def spacing(self):
        """The element spacing. A dict."""
        return self._mesh_._element_spacing_

    @property
    def ratio(self):
        """The element ratio. A dict."""
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
        assert other.__class__.__name__ == '_3dCSCG_Mesh_Elements', \
            'I can only equal to a _3dCSCG_Mesh_Elements object.'
        # when mesh is equal, elements must be
        mesh_judge = self._mesh_ == other._mesh_
        indices_judge = self.indices == other.indices
        judge = mesh_judge and indices_judge
        return cOmm.allreduce(judge, op=MPI.LAND)

    @property
    def quality(self):
        """Return a dict: Keys are element indices, values are the
        qualities of the elements indicated by the keys."""
        if self._quality_ is None:
            self._quality_ = dict()
            tes = self._mesh_.trace.elements
            trace_elements_quality = tes.quality
            for i in self:
                qualities = [0. for _ in range(6)]
                MAP = tes.map[i]
                for j in range(6):
                    qualities[j] = trace_elements_quality[MAP[j]]
                self._quality_[i] = min(qualities)
        return self._quality_

    @property
    def DO(self):
        if self._DO_ is None:
            self._DO_ = _3dCSCG_Mesh_Elements_DO(self)
        return self._DO_

    def ___PRIVATE_elementwise_cache_metric_key___(self, i):
        """A private method to produce str key for cache."""
        mark = self[i].type_wrt_metric.mark
        if not isinstance(mark, str): mark = str(mark)
        return mark

    @accepts('self', int)
    def ___DO_find_slave_of_element___(self, i: int) -> int:
        """Find the core rank of mesh element #i."""
        return self._mesh_.DO.FIND_slave_of_element(i)


class _3dCSCG_Mesh_Elements_DO(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._FIND_ = None
        self._freeze_self_()

    def illustrate_element(self, i, density_factor=2):
        """We use this method to illustrate a mesh element.

        :param int i: We illustrate mesh element #i.
        :param int density_factor: How refined the plots are? be in {1,2,3,...}
        :return:
        """
        if not i in self._elements_:
            cOmm.barrier()
            return
        density = 5 + 4 * density_factor
        i0 = 1 + density_factor
        i1 = 2 * density_factor + 2
        i2 = 3 * density_factor + 3
        _ = np.linspace(-1, 1, density)
        r, s = np.meshgrid(_, _, indexing='ij')
        anchors = ( # the points we will plot the outward unit norm vector
            [i0, i0],
            [i0, i2],
            [i2, i0],
            [i2, i2],
            [i1, i1],
        )
        uv_r = np.array([r[indices[0],indices[1]] for indices in anchors])
        uv_s = np.array([s[indices[0],indices[1]] for indices in anchors])

        element = self._elements_[i]
        sides = element.sides

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        x_lim, y_lim, z_lim = [list(), list(), list()]
        for sn in sides:
            side = sides[sn]
            CT = side.coordinate_transformation
            x, y, z = CT.mapping(r, s)
            ax.plot_surface(x, y, z)
            x_lim.append(np.min(x))
            x_lim.append(np.max(x))
            y_lim.append(np.min(y))
            y_lim.append(np.max(y))
            z_lim.append(np.min(z))
            z_lim.append(np.max(z))

        x_range = np.max(x_lim) - np.min(x_lim)
        y_range = np.max(y_lim) - np.min(y_lim)
        z_range = np.max(z_lim) - np.min(z_lim)

        mean_range = (x_range + y_range + z_range) / 6

        for sn in sides:
            side = sides[sn]
            CT = side.coordinate_transformation
            x, y, z = CT.mapping(uv_r, uv_s)
            u, v, w = CT.outward_unit_normal_vector(uv_r, uv_s)
            ax.quiver(x, y, z, u*mean_range, v*mean_range, w*mean_range, color='red', linewidth=0.8)

        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_zlabel(r'$z$')
        plt.title(f"mesh element: {i}")

        plt.show()
        cOmm.barrier()

    @property
    def FIND(self):
        if self._FIND_ is None:
            self._FIND_ = _3dCSCG_Mesh_Elements_DO_FIND(self._elements_)
        return self._FIND_


class _3dCSCG_Mesh_Elements_DO_FIND(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._freeze_self_()

    def slave_of_element(self, i):
        """Find the core rank of mesh element #i."""
        return self._elements_.___DO_find_slave_of_element___(i)




class _3dCSCG_Mesh_Elements_CT(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self._ctq_1d_ = _3dCSCG_ECT_1d_QUAD(self)
        self._ctq_3d_ = _3dCSCG_ECT_3d_QUAD(self)

    @property
    def QUAD_1d(self):
        """Evaluating the coordinate transformation from quadrature. Results are in 1d array (do ravel)."""
        return self._ctq_1d_

    @property
    def QUAD_3d(self):
        """Evaluating the coordinate transformation from quadrature. Results are in 3d array."""
        return self._ctq_3d_


    def mapping(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta, sigma)
    def X(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta, sigma)
    def Y(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta, sigma)
    def Z(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'Z', xi, eta, sigma)

    def Jacobian_matrix(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta, sigma)
    def J00(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta, sigma)
    def J01(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta, sigma)
    def J02(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J02', xi, eta, sigma)
    def J10(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta, sigma)
    def J11(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta, sigma)
    def J12(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J12', xi, eta, sigma)
    def J20(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J20', xi, eta, sigma)
    def J21(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J21', xi, eta, sigma)
    def J22(self, xi, eta, sigma):
        return ElementsCTValuesCache(self._elements_, 'J22', xi, eta, sigma)

    def Jacobian(self, xi, eta, sigma, J=None):
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, sigma, intermediateData=J)
    def metric(self, xi, eta, sigma, detJ=None):
        """g := det(G) = Jacobian ** 2."""
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, sigma, intermediateData=detJ)
    def metric_matrix(self, xi, eta, sigma, J=None):
        """G, g_{i,j}."""
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, sigma, intermediateData=J)

    def inverse_Jacobian_matrix(self, xi, eta, sigma, J=None):
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, sigma, intermediateData=J)
    def inverse_Jacobian(self, xi, eta, sigma, iJ=None):
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, sigma, intermediateData=iJ)
    def inverse_metric_matrix(self, xi, eta, sigma, iJ=None):
        """G^-1, g^{i,j}."""
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, sigma, intermediateData=iJ)


class ElementsCTValuesCache(FrozenOnly):
    def __init__(self, elements, CTT, xi, eta, sigma, intermediateData=None):
        self._elements_ = elements
        self._CTT_ = CTT
        self._xi_eta_sigma_ = (xi, eta, sigma)
        if intermediateData is not None:
            assert intermediateData.__class__.__name__ == 'ElementsCTValuesCache', \
                "intermediateData can only be an ElementsCTValues object."
            inter_xi_eta_sigma = intermediateData._xi_eta_sigma_
            for i, com in enumerate(inter_xi_eta_sigma):
                assert self._xi_eta_sigma_[i] is com, \
                    "intermediateData xi, eta, sigma must be the same object."
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
        if self._CTT_ in ('mapping', 'X', 'Y', 'Z'):
            return getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
        elif self._CTT_ in ('Jacobian_matrix', 'J00', 'J01', 'J02', 'J10', 'J11', 'J12', 'J20', 'J21', 'J22'):

            if isinstance(type_wrt_metric, int):
                return getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
            elif type_wrt_metric in self._cache_:
                return self._cache_[type_wrt_metric]
            else:
                JM = getattr(element.coordinate_transformation, self._CTT_)(*self._xi_eta_sigma_)
                if type_wrt_metric in self._multi_elements_metric_ and \
                    type_wrt_metric not in self._cache_ and \
                    self._multi_elements_metric_[type_wrt_metric] >= caChe_factor:
                    self._cache_[type_wrt_metric] = JM
                return JM
        else:
            if isinstance(type_wrt_metric, int):
                if self._intermediateData_ is None:
                    return getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_)
                else:
                    return getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_, self._intermediateData_[i])
            elif type_wrt_metric in self._cache_:
                return self._cache_[type_wrt_metric]
            else:
                if self._intermediateData_ is None:
                    result = getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_)
                else:
                    result = getattr(element.coordinate_transformation, self._CTT_)(
                        *self._xi_eta_sigma_, self._intermediateData_[i])
                if type_wrt_metric in self._multi_elements_metric_ and \
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


class _3dCSCG_ECT_1d_QUAD(FrozenOnly):
    def __init__(self, ect):
        self._elements_ = ect._elements_
        self._freeze_self_()

    @staticmethod
    def ___compute_xietasigma___(quad_degree, quad_type):
        _Quadrature_ = Quadrature(quad_degree, category=quad_type)
        quad_nodes = _Quadrature_.quad[0]
        xi, eta, sigma = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        sigma = sigma.ravel('F')
        return xi, eta, sigma

    def mapping(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta, sigma)

    def X(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta, sigma)
    def Y(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta, sigma)
    def Z(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Z', xi, eta, sigma)

    def Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta, sigma)
    def J00(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta, sigma)
    def J01(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta, sigma)
    def J02(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J02', xi, eta, sigma)
    def J10(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta, sigma)
    def J11(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta, sigma)
    def J12(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J12', xi, eta, sigma)
    def J20(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J20', xi, eta, sigma)
    def J21(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J21', xi, eta, sigma)
    def J22(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J22', xi, eta, sigma)

    def Jacobian(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, sigma, intermediateData=None)
    def metric(self, quad_degree, quad_type):
        """g := det(G) = Jacobian ** 2."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, sigma, intermediateData=None)
    def metric_matrix(self, quad_degree, quad_type):
        """G, g_{i,j}."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, sigma, intermediateData=None)

    def inverse_Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, sigma, intermediateData=None)
    def inverse_Jacobian(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, sigma, intermediateData=None)
    def inverse_metric_matrix(self, quad_degree, quad_type):
        """G^-1, g^{i,j}."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, sigma, intermediateData=None)


class _3dCSCG_ECT_3d_QUAD(FrozenOnly):
    def __init__(self, ect):
        self._elements_ = ect._elements_
        self._freeze_self_()

    @staticmethod
    def ___compute_xietasigma___(quad_degree, quad_type):
        _Quadrature_ = Quadrature(quad_degree, category=quad_type)
        quad_nodes = _Quadrature_.quad[0]
        xi, eta, sigma = np.meshgrid(*quad_nodes, indexing='ij')
        return xi, eta, sigma

    def mapping(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'mapping', xi, eta, sigma)

    def X(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'X', xi, eta, sigma)
    def Y(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Y', xi, eta, sigma)
    def Z(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Z', xi, eta, sigma)

    def Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian_matrix', xi, eta, sigma)
    def J00(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J00', xi, eta, sigma)
    def J01(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J01', xi, eta, sigma)
    def J02(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J02', xi, eta, sigma)
    def J10(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J10', xi, eta, sigma)
    def J11(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J11', xi, eta, sigma)
    def J12(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J12', xi, eta, sigma)
    def J20(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J20', xi, eta, sigma)
    def J21(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J21', xi, eta, sigma)
    def J22(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'J22', xi, eta, sigma)

    def Jacobian(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'Jacobian', xi, eta, sigma, intermediateData=None)
    def metric(self, quad_degree, quad_type):
        """g := det(G) = Jacobian ** 2."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric', xi, eta, sigma, intermediateData=None)
    def metric_matrix(self, quad_degree, quad_type):
        """G, g_{i,j}."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'metric_matrix', xi, eta, sigma, intermediateData=None)

    def inverse_Jacobian_matrix(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian_matrix', xi, eta, sigma, intermediateData=None)
    def inverse_Jacobian(self, quad_degree, quad_type):
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_Jacobian', xi, eta, sigma, intermediateData=None)
    def inverse_metric_matrix(self, quad_degree, quad_type):
        """G^-1, g^{i,j}."""
        xi, eta, sigma = self.___compute_xietasigma___(quad_degree, quad_type)
        return ElementsCTValuesCache(self._elements_, 'inverse_metric_matrix', xi, eta, sigma, intermediateData=None)




if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\elements\main.py
    from _3dCSCG.main import MeshGenerator
    elements = [5, 5, 5]
    mesh = MeshGenerator('crazy', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)

    for i in range(mesh.elements.GLOBAL_num):
        mesh.elements.DO.illustrate_element(i)

    # print(mesh.elements.quality)
    # print(mesh.quality)