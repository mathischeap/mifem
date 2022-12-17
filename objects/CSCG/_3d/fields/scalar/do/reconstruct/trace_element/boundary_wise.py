# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class OnTraceElement_BoundaryWise(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel, i):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param i:
            Which trace elements we are going to reconstruct the scalar.

            If `i` is:
                1) one of {None, 'on_mesh_boundaries'}, we reconstruct on all trace-elements on the
                    mesh boundaries.

        :return:
        """
        SELF = self._sf_

        xyz = dict()
        value = dict()

        func = SELF.___DO_evaluate_func_at_time___()

        if i in (None, 'on_mesh_boundaries'):
            RTE = SELF.mesh.boundaries.range_of_trace_elements
            INDICES = list()
            for bn in SELF.func:
                INDICES.extend(RTE[bn])

        else:
            raise NotImplementedError(f"_3dCSCG_ScalarField of ftype 'boundary-wise'"
                                      f"trace-element-reconstruction currently doesn't accept i={i}.")

        for _I in INDICES:
            te = SELF.mesh.trace.elements[_I]
            assert te.IS_on_mesh_boundary, f"must be the case!"
            # when we do trace-element-wise reconstruction, we only accept 1d xi, eta, sigma.
            xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, parse_3_1d_eps=True)

            bn = te.on_mesh_boundary
            assert bn in func, f"trace element #{_I} is on <{bn}> which is not covered by boundary-wise func."
            func_i = func[bn][0]

            v_i = func_i(*xyz_i)

            if ravel:
                xyz[_I] = [_.ravel('F') for _ in xyz_i]
                value[_I] = [v_i.ravel('F'), ]
            else:
                xyz[_I] = xyz_i
                value[_I] = [v_i, ]

        return xyz, value
