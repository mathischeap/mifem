# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly


class _3dCSCG_Mesh_Element_Side_CT(FrozenOnly):
    """The coordinate transformation class for the mesh element side class."""
    def __init__(self, side):
        self._side_ = side
        te = self._side_._mesh_.trace.elements[self._side_.trace_element]
        self._te_ct_ = te.coordinate_transformation
        self._freeze_self_()

    def mapping(self, *evaluation_points):
        return self._te_ct_.mapping(
            *evaluation_points,
            from_element=self._side_._element_.i,
            side=self._side_.side_name,
        )

    def outward_unit_normal_vector(self, *evaluation_points):
        """
        evaluate the outward_unit_normal_vector on points
        *evaluation_points.

        :param evaluation_points:
        :return:
        """
        uv = self._te_ct_.___PRIVATE_outward_unit_normal_vector___(
            *evaluation_points,
            from_element=self._side_._element_.i, side=self._side_.side_name)
        return uv
