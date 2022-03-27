


from screws.warnings.trace_element import TraceElementWarning
from screws.freeze.base import FrozenOnly
import numpy as np
from screws.functions._3d_space.angle import angle_between_two_vectors

import warnings


class _3dCSCG_Trace_Elements_SELFCHECK(FrozenOnly):
    """We find some specific groups of elements."""

    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._freeze_self_()

    def outward_unit_normal_vector(self):
        """This is a self check program to check that we obtain the
        correct outward unit normal vector(s) for a trace element.

        Notice this self check program does not give error. It will only
        give warning; it can only find something may be wrong.

        We basically compute the angle between the outward norm vector
        and the vector pointing the mesh element center. If the angle
        is pi (180 degree), of course, the outward normal vector is
        correct. If it is 0, then we have a wrong outward normal vector.
        """
        # the center of the mesh element.
        SELF = self._elements_
        xi, et, sg = np.array([0,]), np.array([0,]), np.array([0,])
        # the center of the trace element.
        rc, sc = np.array([0,]), np.array([0,])
        for i in SELF:
            te = SELF[i]
            e1 = te.CHARACTERISTIC_element
            s1 = te.CHARACTERISTIC_side
            E = [e1,]
            S = [s1,]
            if not te.IS_on_mesh_boundary:
                p2 = te.NON_CHARACTERISTIC_position
                e2 = int(p2[:-1])
                s2 = p2[-1]
                if e2 in SELF._mesh_.elements:
                    E.append(e2)
                    S.append(s2)
            for e, s in zip(E, S):
                me = SELF._mesh_.elements[e]
                me_ct = me.coordinate_transformation.mapping(xi, et, sg)
                te_ct = te.coordinate_transformation.mapping(rc, sc, from_element=e, side=s)
                te_uv = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(rc, sc, from_element=e, side=s)

                v_inner = (me_ct[0]-te_ct[0], me_ct[1]-te_ct[1], me_ct[2]-te_ct[2])
                v_outer = te_uv
                angle = angle_between_two_vectors(v_inner, v_outer)

                if angle < np.pi/4:
                    warnings.warn(
                        f"___PRIVATE_outward_unit_normal_vector___ of trace element "
                        f"#{i} may be wrong",
                        TraceElementWarning)
                elif angle < np.pi/2:
                    warnings.warn(
                        f"the mesh element #{e} may be very distorted.",
                        TraceElementWarning)
                else:
                    pass
