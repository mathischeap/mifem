

from screws.frozen import FrozenOnly
import numpy as np
from _2dCSCG.mesh.elements.element.coordinate_transformation import _2dCSCG_Mesh_ECT



class _2dCSCG_Mesh_Element(FrozenOnly):
    """"""
    def __init__(self, elements, i):
        self._elements_ = elements
        self._mesh_ = elements._mesh_
        self._i_ = i
        self._type_wrt_metric_ = None
        self._in_region_ = self._mesh_.do.find.region_name_of_element(self.i)
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
        region, localRegionIndices = self._mesh_.do.find.region_name_and_local_indices_of_element(self.i)
        elementsSpacing = self._elements_.spacing[region]
        _spacing_ = np.zeros((2,2))
        for i in range(2):
            _spacing_[i, 0] = elementsSpacing[i][localRegionIndices[i]]
            _spacing_[i, 1] = elementsSpacing[i][localRegionIndices[i]+1]
        return _spacing_

    @property
    def type_wrt_metric(self):
        if self._type_wrt_metric_ is None:
            region, _ = self._mesh_.do.find.region_name_and_local_indices_of_element(self.i)
            self._type_wrt_metric_ = \
                self._mesh_.domain.regions[region].type_wrt_metric.___CLASSIFY_ELEMENT_of_spacing___(
                    self.spacing)
        return self._type_wrt_metric_

    @property
    def coordinate_transformation(self):
        return self._ct_


