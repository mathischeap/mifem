# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('/')

from screws.freeze.main import FrozenOnly
from root.config.main import *
from objects.CSCG._3d.mesh.elements.element.main import _3dCSCG_Mesh_Element
from objects.CSCG._3d.mesh.elements.coordinate_transformation.main import _3dCSCG_Mesh_Elements_CT
from objects.CSCG._3d.mesh.elements.do.main import _3dCSCG_Mesh_Elements_DO
from objects.CSCG._3d.mesh.elements.visualize import _3dCSCG_MeshElements_VIS
from objects.CSCG._3d.mesh.elements.IS import _3dCSCG_MeshElements_IS
from objects.CSCG._3d.mesh.elements.find import _3dCSCG_MeshElements_Find


class _3dCSCG_Mesh_Elements(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = dict()
        self._DO_ = None
        self._visualize_ = None
        self._IS_ = None
        self._find_ = None
        self._ct_ = _3dCSCG_Mesh_Elements_CT(self)
        for i in self.indices:
            self._elements_[i] = _3dCSCG_Mesh_Element(self, i)
        self.___PRIVATE_parse_elements_type_wrt_metric___()
        self.___PRIVATE_reset_cache___()
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
        """
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
        Statistic['similarity according to types wrt metric'] = similarity  # in [0, 1]
        Statistic['num_local_orthogonal_elements'] = self._num_local_orthogonal_elements_

        return Statistic

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

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _3dCSCG_MeshElements_Find(self)
        return self._find_






if __name__ == '__main__':
    # mpiexec -n 4 python _3dCSCG\mesh\elements\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [1, 1, 1]
    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1], [0,1], [0,1]))(elements)

    # for i in range(mesh.elements.GLOBAL_num):
    #     mesh.elements.do.illustrate_element(i)

    print(mesh.elements.statistic)
    # print(mesh.quality)