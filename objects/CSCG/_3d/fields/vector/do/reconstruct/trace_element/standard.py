# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly


class OnTraceElement_Standard(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel, i):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param i:
            1) self.ftype == 'standard':
                Do the reconstruction in mesh element #i. When it is None, it means all local mesh
                elements.
        :return:
        """

        SELF = self._vf_

        xyz = dict()
        value = dict()

        func = SELF.___DO_evaluate_func_at_time___()


        if isinstance(i, int):
            INDICES = [i,]
        elif i == 'on_mesh_boundaries': # then we plot on all mesh boundaries (mesh elements on the boundaries)
            INDICES = list()
            RTE = SELF.mesh.boundaries.range_of_trace_elements
            for bn in RTE:
                INDICES.extend(RTE[bn])
        else:
            raise NotImplementedError(f"_3dCSCG_VectorField of 'standard' ftype"
                                      f" trace-element-reconstruction currently don't accept i={i}.")

        for I in INDICES:
            te = SELF.mesh.trace.elements[I]
            # when reconstruct on trace-element-wise, we only accept 1d xi, eta and sigma.
            xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, parse_3_1d_eps=True)

            vx_i = func[0](*xyz_i)
            vy_i = func[1](*xyz_i)
            vz_i = func[2](*xyz_i)

            if ravel:
                xyz[I] = [_.ravel('F') for _ in xyz_i]
                value[I] = [vx_i.ravel('F'), vy_i.ravel('F'), vz_i.ravel('F')]
            else:
                xyz[I] = xyz_i
                value[I] = [vx_i, vy_i, vz_i,]

        return xyz, value