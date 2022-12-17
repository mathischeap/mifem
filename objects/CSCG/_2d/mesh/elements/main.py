# -*- coding: utf-8 -*-
import sys
if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly

from objects.CSCG._2d.mesh.elements.coordinate_transformation.main import _2dCSCG_Mesh_Elements_CT
from root.config.main import *
from objects.CSCG._2d.mesh.elements.element.main import _2dCSCG_Mesh_Element
from objects.CSCG._2d.mesh.elements.whether import _2dCSCG_MeshElements_Whether
from objects.CSCG._2d.mesh.elements.visualize import _2dCSCG_MeshElements_VIS
from objects.CSCG._2d.mesh.elements.do.main import _2dCSCG_Mesh_Elements_do


class _2dCSCG_Mesh_Elements(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = dict()
        self._whether_ = None
        self._visualize_ = None
        self._do_ = None
        self._ct_ = _2dCSCG_Mesh_Elements_CT(self)
        for i in self.indices:
            self._elements_[i] = _2dCSCG_Mesh_Element(self, i)
        self.___PRIVATE_parse_elements_type_wrt_metric___()
        self._rwPc_ = None
        self._freeze_self_()

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
        if self.num == 0:
            similarity = 'nan'
        elif self.num == 1:
            similarity = 1
        else:
            similarity = (self.num - A) / (self.num - 1)
        Statistic['similarity according to types wrt metric'] = similarity # in [0, 1]
        Statistic['num_local_orthogonal_elements'] = self._num_local_orthogonal_elements_

        return Statistic

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = _2dCSCG_MeshElements_Whether(self)
        return self._whether_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2dCSCG_MeshElements_VIS(self)
        return self._visualize_

    @property
    def global_num(self):
        """The total amount of mesh elements in all cores."""
        return self._mesh_._num_total_elements_

    @property
    def in_regions(self):
        """The regions the local mesh elements of this core are in."""
        return self._mesh_._elements_in_regions_

    @property
    def region_wise_prime_core(self):
        """A dict same in all cores. It shows which core (value) has the most elements of
        a particular region (key).

        If two cores have the same amount of mesh-elements of a particular region, we use the
        core of lower rank as the prime core.

        """
        if self._rwPc_ is None:

            rns = self._mesh_.domain.regions.names
            COUNT = dict()
            for rn in rns:
                COUNT[rn] = 0

            for i in self:
                rn = self._mesh_.do.find.region_name_of_element(i)
                COUNT[rn] += 1

            CTs = COMM.gather(COUNT, root=MASTER_RANK)

            if RANK == MASTER_RANK:
                overall_count = dict()
                for rn in rns:
                    overall_count[rn] = [0 for _ in range(SIZE)]

                for i, CT in enumerate(CTs):
                    for rn in CT:
                        overall_count[rn][i] = CT[rn]

                self._rwPc_ = dict()
                for rn in overall_count:
                    Pc = overall_count[rn].index(max(overall_count[rn]))
                    self._rwPc_[rn] = Pc

            else:
                pass

            self._rwPc_ = COMM.bcast(self._rwPc_, root=MASTER_RANK)

        return self._rwPc_

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
        return COMM.allreduce(judge, op=MPI.LAND)

    def ___PRIVATE_elementwise_cache_metric_key___(self, i):
        """This is an import method, for example, it can be used as a key generator for EWC objects."""
        mark = self[i].type_wrt_metric.mark
        if not isinstance(mark, str): mark = str(mark)
        return mark

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _2dCSCG_Mesh_Elements_do(self)
        return self._do_


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/mesh/elements/main.py
    from objects.CSCG._2d.master import MeshGenerator

    # mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    mesh = MeshGenerator('crazy', c=0, bounds=([0, 1],[0, 1]))([10, 10])
    elements = mesh.elements
    S = elements.statistic
    # print(S)
    print(elements.whether.homogeneous_according_to_types_wrt_metric)
