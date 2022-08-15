# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from root.config.main import *
from objects.CSCG._3d.mesh.elements.element.main import _3dCSCG_Mesh_Element
from objects.CSCG._3d.mesh.elements.coordinate_transformation.main import _3dCSCG_Mesh_Elements_CT
from objects.CSCG._3d.mesh.elements.do.main import _3dCSCG_Mesh_Elements_DO
from objects.CSCG._3d.mesh.elements.visualize import _3dCSCG_MeshElements_VIS
from objects.CSCG._3d.mesh.elements.IS import _3dCSCG_MeshElements_IS


class _3dCSCG_Mesh_Elements(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = dict()
        self._DO_ = None
        self._visualize_ = None
        self._IS_ = None
        self._ct_ = _3dCSCG_Mesh_Elements_CT(self)
        for i in self.indices:
            self._elements_[i] = _3dCSCG_Mesh_Element(self, i)
        self.___PRIVATE_parse_elements_type_wrt_metric___()
        self.___PRIVATE_reset_cache___()
        self._Statistic_ = None
        self._involved_mesh_boundaries_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self.coordinate_transformation.___PRIVATE_reset_cache___()
        self._quality_ = None

    def ___PRIVATE_parse_elements_type_wrt_metric___(self):
        counter: dict = dict()
        self._multi_elements_metric_: dict = dict()
        self._num_local_orthogonal_elements_ = 0
        for i in self:
            ei = self[i]
            mki = ei.type_wrt_metric.mark

            if isinstance(mki, str) and mki[:4] == 'Orth':
                self._num_local_orthogonal_elements_ += 1

            if mki in counter:
                ei._type_wrt_metric_ = self[counter[mki]]._type_wrt_metric_
                if mki in self._multi_elements_metric_:
                    self._multi_elements_metric_[mki] += 1
                else:
                    self._multi_elements_metric_[mki] = 2
            else:
                counter[mki] = i

    @property
    def statistic(self):
        """The statistic of local mesh elements.

        :return: a dict of the statistic.
            Keys:
                - 'total num of local mesh elements':
                - 'amount of types wrt metric':
                - 'similarity according to types wrt metric': [0,1]; if there is no local mesh elements,
                    it is 'nan'
                - 'num_local_orthogonal_elements':
                - 'amount of internal mesh elements': amount of mesh elements none of whose 6 sides
                    is on the mesh boundary. So, it is possible that it has edge on the mesh
                    boundary. For example, consider an L-shape domain, the corner mesh-element
                    will have an edge on the corner of L-shape domain. Thus, thus this edge
                    (edge-element in fact) is on the mesh boundary.
                - 'amount of boundary mesh elements': amount of mesh elements who at least have one
                    side (trace-element) being on the mesh-boundary.
        """
        if self._Statistic_ is None:

            MEM = self._multi_elements_metric_
            Statistic = dict()
            Statistic['total num of local mesh elements'] = self.num
            A = self.num
            for type_name in MEM:
                A -= MEM[type_name] - 1
            Statistic['amount of types wrt metric'] = A
            if self.num == 0:
                similarity = 'nan'
            elif self.num == 1:
                similarity = 1
            else:
                similarity = (self.num - A) / (self.num - 1)
            Statistic['similarity according to types wrt metric'] = similarity  # in [0, 1] or 'nan'
            Statistic['num_local_orthogonal_elements'] = self._num_local_orthogonal_elements_

            AoIME = 0
            AoBME = 0
            for i in self:
                element = self[i]
                if element.IS.internal:
                    AoIME += 1
                else:
                    AoBME += 1
            # note that internal mesh element could have edges on the mesh boundary.
            Statistic['amount of internal mesh elements'] = AoIME
            Statistic['amount of boundary mesh elements'] = AoBME

            self._Statistic_ = Statistic
        return self._Statistic_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _3dCSCG_MeshElements_IS(self)
        return self._IS_

    @property
    def GLOBAL_num(self):
        """The total number of elements in all cores."""
        return self._mesh_._num_total_elements_

    @property
    def in_regions(self):
        """The elements in this core are in these regions."""
        return self._mesh_._elements_in_regions_

    @property
    def involved_mesh_boundaries(self):
        """{List[str]}: Return a list of mesh boundary names on which local mesh-elements has side.
        """
        if self._involved_mesh_boundaries_ is None:
            RoES = mesh.boundaries.range_of_element_sides
            _involved_mesh_boundaries_ = list()
            for bn in RoES:
                if len(RoES[bn]) > 0:
                    _involved_mesh_boundaries_.append(bn)
                else:
                    pass
            self._involved_mesh_boundaries_ = _involved_mesh_boundaries_
        return self._involved_mesh_boundaries_

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
        """return the local element #i"""
        return self._elements_[i]

    def __contains__(self, i):
        """If element #i is a local element?"""
        return i in self.indices

    def __iter__(self):
        """Go through all local elements"""
        for i in self.indices:
            yield i

    def __len__(self):
        """How many local elements?"""
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
        """Return a dict: Keys are element indices, values are the qualities of the elements
        indicated by the keys.
        """
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
    def do(self):
        if self._DO_ is None:
            self._DO_ = _3dCSCG_Mesh_Elements_DO(self)
        return self._DO_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_MeshElements_VIS(self)
        return self._visualize_

    def ___PRIVATE_elementwise_cache_metric_key___(self, i):
        """A private method to produce str key for cache."""
        mark = self[i].type_wrt_metric.mark
        if not isinstance(mark, str): mark = str(mark)
        return mark








if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\mesh\elements\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [1, 1, 1]
    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1], [0,1], [0,1]))(elements)

    # for i in range(mesh.elements.GLOBAL_num):
    #     mesh.elements.do.illustrate_element(i)

    print(mesh.elements.statistic)
    # print(mesh.quality)