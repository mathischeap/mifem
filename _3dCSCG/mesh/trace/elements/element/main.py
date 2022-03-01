

import sys
if './' not in sys.path: sys.path.append('./')

from root.config import rAnk


from screws.frozen import FrozenOnly

from _3dCSCG.mesh.trace.elements.element.coordinate_transformation import _3dCSCG_Trace_Element_CoordinateTransformation




class _3dCSCG_Trace_Element(FrozenOnly):
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
            "CHARACTERISTIC_element must be in the same core."
        if self.IS_on_mesh_boundary:
            assert self.NON_CHARACTERISTIC_position[0] not in '1234567890'
        self._ct_ = None
        self._type_wrt_metric_ = None
        self._freeze_self_()
        # # do a check for periodic trace element ________________________
        # if self.IS_on_periodic_boundary:
        #     assert not self.IS_on_mesh_boundary # must be the case
        #     e1 = int(position_1[:-1])
        #     e2 = int(position_2[:-1])
        #     if e1 == e2:
        #         warnings.warn(f"periodic trace element #{i} is on two "
        #                       f"faces of same mesh element #{e1}, "
        #                       f"this may cause unknown error. Please "
        #                       f"consider use more mesh elements to "
        #                       f"remove this situation.",
        #                       PeriodicTraceElementWarning)

    @property
    def positions(self):
        return self._p1_, self._p2_

    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = _3dCSCG_Trace_Element_CoordinateTransformation(self)
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
        """The other position."""
        return self._ncp_

    @property
    def CHARACTERISTIC_position(self):
        """The position we mainly locate this trace element."""
        return self._cp_
    @property
    def CHARACTERISTIC_element(self):
        """We mainly consider this trace element is a side of this mesh
        element."""
        return int(self._cp_[:-1])
    @property
    def CHARACTERISTIC_side(self):
        """We main consider this trace element is such a side of the
        CHARACTERISTIC_element."""
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
        side = self.CHARACTERISTIC_side
        if side == 'N':
            trace_spacing = (element_spacing[0][0], element_spacing[1], element_spacing[2])
        elif side == 'S':
            trace_spacing = (element_spacing[0][1], element_spacing[1], element_spacing[2])
        elif side == 'W':
            trace_spacing = (element_spacing[0], element_spacing[1][0], element_spacing[2])
        elif side == 'E':
            trace_spacing = (element_spacing[0], element_spacing[1][1], element_spacing[2])
        elif side == 'B':
            trace_spacing = (element_spacing[0], element_spacing[1], element_spacing[2][0])
        elif side == 'F':
            trace_spacing = (element_spacing[0], element_spacing[1], element_spacing[2][1])
        else:
            raise Exception()
        return trace_spacing

    @property
    def i(self):
        """This is the ith trace element."""
        return self._i_

    @property
    def IS_on_mesh_boundary(self):
        """As this property name says."""
        return self._ondb_

    @property
    def on_mesh_boundary(self):
        """Return the mesh boundary name this trace element is on. If it is not on one, return None."""
        if self.IS_on_mesh_boundary:
            return self.NON_CHARACTERISTIC_position
        else:
            return None

    @property
    def IS_on_periodic_boundary(self):
        """As this property name says."""
        return self._onpb_

    @property
    def IS_shared_by_cores(self):
        """True or False, as this property name says."""
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
        """None or a int. """
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

    @property
    def type_wrt_metric(self):
        """Return the trace-element-metric-type object reflecting the element type."""
        if self._type_wrt_metric_ is None:

            self._type_wrt_metric_ = \
                self._mesh_.domain.regions[
                    self.CHARACTERISTIC_region].type_wrt_metric.___CLASSIFY_TRACE_ELEMENT_of_spacing___(
                    self.spacing)

        return self._type_wrt_metric_







if __name__ == '__main__':
    # mpiexec -n 12 python _3dCSCG\mesh\trace\elements\element\main.py
    from _3dCSCG.main import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.SELFCHECK.outward_unit_normal_vector()
    Q = mesh.trace.elements.quality
    print(mesh.quality)
    print(mesh.trace.quality)

    mesh.trace.elements.do.illustrate_trace_element(1)

    # te0 = mesh.trace.elements[0]

    # print(te0.IS_on_periodic_boundary)

    for i in range(mesh.trace.elements.GLOBAL_num):
        if i in mesh.trace.elements:
            te = mesh.trace.elements[i]

            print(rAnk, te.type_wrt_metric.mark)