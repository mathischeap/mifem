# -*- coding: utf-8 -*-
""""""
import sys
if './' not in sys.path: sys.path.append('./')

from components.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.elements.element.sub_geometry.sub_geometry import ElementSubGeometry
import numpy as np
from objects.CSCG._3d.mesh.elements.element.sides.main import _3dCSCG_Mesh_Element_Sides
from objects.CSCG._3d.mesh.elements.element.coordinate_transformation import _3dCSCG_Mesh_Element_CT
from objects.CSCG._3d.mesh.elements.element.do import _3dCSCG_MeshElement_Do
from objects.CSCG._3d.mesh.elements.element.IS import _3dCSCG_MeshElement_IS

class _3dCSCG_Mesh_Element(FrozenOnly):
    """The mesh element class"""
    def __init__(self, elements, i):
        self._elements_ = elements
        self._mesh_ = elements._mesh_
        self._i_ = i
        self._type_wrt_metric_ = None
        self._in_region_ = self._mesh_.do.find.region_name_of_element(i)
        self._ct_ = None
        self._sub_geometry_ = None
        self._sides_ = None
        self._do_ = None
        self._IS_ = None
        self._freeze_self_()

    @property
    def i(self):
        """The global numbering of this element."""
        return self._i_

    @property
    def position(self):
        """The elements.map[i] reflects the position of an element."""
        return self._elements_.map[self.i]

    @property
    def in_region(self):
        """This element is in which domain regions?"""
        return self._in_region_

    @property
    def region_indices(self):
        """This element's indices in the region is this one."""
        return self._mesh_.do.find.region_name_and_local_indices_of_element(self.i)[1]

    @property
    def spacing(self):
        """What is the spacing of this element in the domain regions?

        This property basically reflects the relative position of this element in the domain regions.
        """
        region, localRegionIndices = self._mesh_.do.find.region_name_and_local_indices_of_element(self.i)
        elementsSpacing = self._elements_.spacing[region]
        _spacing_ = np.zeros((3,2))
        for i in range(3):
            _spacing_[i, 0] = elementsSpacing[i][localRegionIndices[i]]
            _spacing_[i, 1] = elementsSpacing[i][localRegionIndices[i]+1]
        return _spacing_

    @property
    def type_wrt_metric(self):
        """Return an element metric type object reflecting the element type."""
        if self._type_wrt_metric_ is None:
            self._type_wrt_metric_ = \
                self._mesh_.domain.regions[
                    self.in_region].type_wrt_metric.___CLASSIFY_ELEMENT_of_spacing___(
                    self.spacing)
        return self._type_wrt_metric_

    @property
    def coordinate_transformation(self):
        """The coordinate transformation object of this element."""
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Mesh_Element_CT(self)
        return self._ct_

    @property
    def sub_geometry(self):
        if self._sub_geometry_ is None:
            self._sub_geometry_ = ElementSubGeometry(self)
        return self._sub_geometry_

    @property
    def sides(self):
        if self._sides_ is None:
            self._sides_ = _3dCSCG_Mesh_Element_Sides(self)
        return self._sides_

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _3dCSCG_MeshElement_Do(self)
        return self._do_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _3dCSCG_MeshElement_IS(self)
        return self._IS_








if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\mesh\elements\element\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [2, 2, 2]
    mesh = MeshGenerator('crazy', c=0.3, bounds=([0,3], [0,3], [0,3]))(elements)

    if 0 in mesh.elements:
        e = mesh.elements[0]
        ess = e.sides
        N = ess['S']
        print(N.trace_element)