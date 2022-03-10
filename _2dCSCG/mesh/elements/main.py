# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly

from _2dCSCG.mesh.elements.coordinate_transformation import _2dCSCG_Mesh_Elements_CT
from root.config.main import *
from _2dCSCG.mesh.elements.element.main import _2dCSCG_Mesh_Element


class _2dCSCG_Mesh_Elements(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = dict()
        self._ct_ = _2dCSCG_Mesh_Elements_CT(self)
        for i in self.indices:
            self._elements_[i] = _2dCSCG_Mesh_Element(self, i)
        self.___PRIVATE_parse_elements_type_wrt_metric___()
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        """"""
        self.coordinate_transformation.___PRIVATE_reset_cache___()

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
    def statistic(self):
        """The statistic of local mesh elements according to their `type_wrt_metric`.

        :return: a dict of the statistic
        """
        MEM = self._multi_elements_metric_
        Statistic = dict()
        Statistic['total num of local mesh elements'] = self.num
        A = self.num
        for type_name in MEM:
            A -= MEM[type_name] - 1
        Statistic['amount of types wrt metric'] = A
        Statistic['similarity according to types wrt metric'] = (self.num-A)/(self.num-1) # in [0, 1]

        return Statistic


    @property
    def GLOBAL_num(self):
        """The total amount of mesh elements in all cores."""
        return self._mesh_._num_total_elements_

    @property
    def in_regions(self):
        """The regions the local mesh elements of this core are in."""
        return self._mesh_._elements_in_regions_

    @property
    def indices(self):
        """The indices (numbers) of all local mesh elements."""
        return self._mesh_._element_indices_

    @property
    def num(self):
        """The amount of local mesh elements."""
        return self._mesh_._num_local_elements_

    @property
    def map(self):
        """The element map of local mesh elements."""
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
        """This is an import method, for example, it can be used as a key generator for EWC objects."""
        mark = self[i].type_wrt_metric.mark
        if not isinstance(mark, str): mark = str(mark)
        return mark





if __name__ == '__main__':
    # mpiexec python _2dCSCG\mesh\elements\main.py
    from _2dCSCG.main import MeshGenerator

    # mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    mesh = MeshGenerator('crazy', c=0.1, bounds=([0,1],[0,1]))([10,10])
    elements = mesh.elements
    S = elements.statistic
    print(S)
