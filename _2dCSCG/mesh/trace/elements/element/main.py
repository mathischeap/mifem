

from screws.freeze.main import FrozenOnly

from _2dCSCG.mesh.trace.elements.element.coordinate_transformation import _2dCSCG_Trace_Element_CoordinateTransformation
from root.config.main import rAnk



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
        if self.IS_on_mesh_boundary:
            assert self.NON_CHARACTERISTIC_position[0] not in '1234567890'
        self._ct_ = _2dCSCG_Trace_Element_CoordinateTransformation(self)
        self._freeze_self_()

    @property
    def positions(self):
        return self._p1_, self._p2_

    @property
    def coordinate_transformation(self):
        return self._ct_

    @property
    def normal_direction(self):
        """"""
        if self._p1_[-1] in 'NS':
            return 'NS'
        elif self._p1_[-1] in 'WE':
            return 'WE'
        elif self._p1_[-1] in 'BF':
            return 'BF'
        else:
            raise Exception()

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
    def i(self):
        return self._i_

    @property
    def IS_on_mesh_boundary(self):
        return self._ondb_

    @property
    def IS_on_periodic_boundary(self):
        return self._onpb_

    @property
    def IS_shared_by_cores(self):
        """"""
        if self.IS_on_mesh_boundary:
            return False
        else:
            if int(self._p1_[:-1]) in self._elements_._mesh_.elements and \
                int(self._p2_[:-1]) in self._elements_._mesh_.elements:
                return False
            else:
                return True

    @property
    def shared_with_core(self):
        if self.IS_shared_by_cores:
            if int(self._p1_[:-1]) in self._elements_._mesh_.elements:
                CORE = self._elements_._mesh_.do.find.slave_of_element(int(self._p2_[:-1]))
            elif int(self._p2_[:-1]) in self._elements_._mesh_.elements:
                CORE = self._elements_._mesh_.do.find.slave_of_element(int(self._p1_[:-1]))
            else:
                raise Exception()
            assert CORE != rAnk
            return CORE
        else:
            return None