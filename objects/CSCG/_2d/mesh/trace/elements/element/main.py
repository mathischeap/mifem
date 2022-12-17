# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly

from objects.CSCG._2d.mesh.trace.elements.element.whether import _2dCSCG_TraceElement_whether
from objects.CSCG._2d.mesh.trace.elements.element.coordinate_transformation.main import \
    _2dCSCG_Trace_Element_CoordinateTransformation
from root.config.main import RANK


class _2dCSCG_Trace_Element(FrozenOnly):
    """

    :param trace_elements:
    :param i:
    :param position_1:
    :param position_2:
    :param cp: characteristic position
    :param ondb:
    :param onpb:
    """
    def __init__(self, trace_elements, i, position_1, position_2, cp, ondb=False, onpb=False):
        self._elements_ = trace_elements
        self._mesh_ = trace_elements._mesh_
        self._i_ = i
        self._p1_ = position_1
        self._p2_ = position_2
        self._cp_ = cp
        if position_1 == cp:
            self._ncp_ = position_2
        elif position_2 == cp:
            self._ncp_ = position_1
        else:
            raise Exception()
        self._ondb_ = ondb
        self._onpb_ = onpb
        assert self.CHARACTERISTIC_element in self._elements_._mesh_.elements, \
            "CHARACTERISTIC_element must be int the same core."
        if self._ondb_:
            assert self.NON_CHARACTERISTIC_position[0] not in '1234567890'
        self._ct_ = None
        self._whether_ = None
        self._type_wrt_metric_ = None
        self._freeze_self_()

    @property
    def positions(self):
        return self._p1_, self._p2_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _2dCSCG_Trace_Element_CoordinateTransformation(self)
        return self._ct_

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = _2dCSCG_TraceElement_whether(self)
        return self._whether_

    @property
    def normal_direction(self):
        """"""
        if self._p1_[-1] in 'UD':
            return 'UD'
        elif self._p1_[-1] in 'LR':
            return 'LR'
        else:
            raise Exception(self._p1_)

    @property
    def on_mesh_boundary(self):
        """Return the mesh boundary name this trace element is on. If it is not on one, return None."""
        if self._ondb_:
            return self.NON_CHARACTERISTIC_position
        else:
            return None

    @property
    def NON_CHARACTERISTIC_position(self):
        return self._ncp_

    @property
    def CHARACTERISTIC_position(self):
        return self._cp_

    @property
    def CHARACTERISTIC_element(self):
        return int(self._cp_[:-1])

    @property
    def CHARACTERISTIC_edge(self):
        return self._cp_[-1]

    @property
    def CHARACTERISTIC_region(self):
        """We mainly consider this trace element is in this regions."""
        region = self._mesh_.do.find.region_name_of_element(
            self.CHARACTERISTIC_element)
        return region

    @property
    def spacing(self):
        element_spacing = self._mesh_.elements[self.CHARACTERISTIC_element].spacing
        edge = self.CHARACTERISTIC_edge
        if edge == 'U':
            trace_spacing = (element_spacing[0][0], element_spacing[1])
        elif edge == 'D':
            trace_spacing = (element_spacing[0][1], element_spacing[1])
        elif edge == 'L':
            trace_spacing = (element_spacing[0], element_spacing[1][0])
        elif edge == 'R':
            trace_spacing = (element_spacing[0], element_spacing[1][1])
        else:
            raise Exception()
        return trace_spacing

    @property
    def i(self):
        return self._i_

    @property
    def shared_with_core(self):
        if self.whether.shared_by_cores:
            if int(self._p1_[:-1]) in self._elements_._mesh_.elements:
                CORE = self._elements_._mesh_.do.find.slave_of_element(int(self._p2_[:-1]))
            elif int(self._p2_[:-1]) in self._elements_._mesh_.elements:
                CORE = self._elements_._mesh_.do.find.slave_of_element(int(self._p1_[:-1]))
            else:
                raise Exception()
            assert CORE != RANK
            return CORE
        else:
            return None

    @property
    def type_wrt_metric(self):
        """Return the trace-element-metric-type object reflecting the element type."""
        if self._type_wrt_metric_ is None:

            self._type_wrt_metric_ = \
                self._mesh_.domain.regions[
                    self.CHARACTERISTIC_region].type_wrt_metric.___CLASSIFY_TRACE_ELEMENT_of_spacing___(
                    self.spacing)

        return self._type_wrt_metric_
