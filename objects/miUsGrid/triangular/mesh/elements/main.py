# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/16 9:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.mesh.elements.do.main import miUsGrid_TriangularMesh_Elements_DO
from objects.miUsGrid.triangular.mesh.elements.num import miUsGrid_TriangularMesh_Elements_Num
from objects.miUsGrid.triangular.mesh.elements.element.main import miUsGrid_TriangularMesh_Element
from objects.miUsGrid.triangular.mesh.elements.coordinate_transformation.main import miUsTriangle_Elements_CoordinateTransformation


class miUsGrid_TriangularMesh_Elements(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._range_ = None
        self._distributions_ = None
        self._map_ = None
        self._cells_ = None
        self._points_ = None
        self.__num_GLOBAL_points__ = None
        self._do_ = miUsGrid_TriangularMesh_Elements_DO(self)
        self._num_ = miUsGrid_TriangularMesh_Elements_Num(self)
        self._elements_ = dict()
        self._similarity_ = None
        self.___hSkeyG___ = None
        self._ct_ = miUsTriangle_Elements_CoordinateTransformation(self)
        self._freeze_self_()

    @property
    def map(self):
        """The element map of all local cells."""
        return self._map_

    @property
    def distributions(self):
        """A list of shape `(sIze,)`.

        How many local cells in each core? distribution[i] means core #`i` have this many local
         cells.
         """
        return self._distributions_

    @property
    def range(self):
        """The local cells numbers are exactly elements of this range."""
        return self._range_

    @property
    def indices(self):
        """Same as range."""
        return self._range_

    @property
    def cells(self):
        """A dict whose keys are local cells, values are corresponding vertexes (points)."""
        return self._cells_

    @property
    def points(self):
        """A dict, keys are numbers of local points (all the vertexes),
        values are corresponding coordinates."""
        return self._points_

    def __iter__(self):
        """Go through all local cells."""
        for i in self.range:
            yield i

    def __contains__(self, item):
        """If item is a local cell number?"""
        return item in self.range

    def __len__(self):
        """How many local cells?"""
        return len(self.range)

    def __getitem__(self, cell_number):
        """"""
        if cell_number not in self._elements_:
            self._elements_[cell_number] = miUsGrid_TriangularMesh_Element(self, cell_number)
        return self._elements_[cell_number]

    @property
    def do(self):
        return self._do_

    @property
    def num(self):
        return self._num_

    @property
    def similarity(self):
        return self._similarity_

    def ___Pr_analyze_element_shapes___(self):
        """We will find shape indicators of all elements and check if it's worthy to store those
        shape indicators for caching. For example, if almost all elements are of different shapes,
        then we just directly consider all elements to be different. If these are quite a big amount of
        elements are of same shapes, we then could use the shape indicators as keys to do the cache.
        """
        shape_indicator_dict = set()
        num_types = 0
        num_elements = self.num.cells
        for e in self:
            element = self[e]
            esi = element.shape_indicator
            if esi not in shape_indicator_dict:
                num_types += 1
                shape_indicator_dict.add(esi)
        if num_elements <= 1:
            similarity = 0
        else:
            factor = num_types / num_elements
            MIN = 1 / num_elements
            MAX = 1

            similarity = (1 - factor) / (MAX - MIN)
        self._similarity_ = similarity

    @property
    def ___Pr_EWC_cache_key___(self):
        if self._similarity_ >= 0.5:
            return self.___Pr_EWC_high_similarity_key_generator___
        else:
            return 'no_cache'

    def ___Pr_EWC_high_similarity_key_generator___(self, i):
        """"""
        if self.___hSkeyG___ is None:
            kG = dict()
            for e in self:
                element = self[e]
                esi = element.shape_indicator
                if esi not in kG:
                    kG[esi] = [e,]
                else:
                    kG[esi].append(e)

            self.___hSkeyG___ = dict()
            for esi in kG:
                for e in kG[esi]:
                    self.___hSkeyG___[e] = esi

        return self.___hSkeyG___[i]

    @property
    def coordinate_transformation(self):
        return self._ct_




if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/elements/main.py

    from __init__ import miTri
    fc = miTri.form('rand0', 2)
    mesh = fc.mesh
    print(mesh.elements.similarity)

    # for e in mesh.elements:
    #     element = mesh.elements[e]
