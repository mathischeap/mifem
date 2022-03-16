from screws.freeze.inheriting.frozen_only import FrozenOnly


class OnTraceElement_BoundaryWise(FrozenOnly):
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

        if i in (None, 'on_mesh_boundaries'):
            INDICES = list()
            RTE = SELF.mesh.boundaries.range_of_trace_elements
            for bn in SELF.func:  # this may not contain all mesh boundaries, only valid ones.
                INDICES.extend(RTE[bn])

        else:
            raise NotImplementedError(f"_3dCSCG_VectorField of 'boundary-wise' ftype"
                                      f" trace-element-reconstruction currently don't accept i={i}.")

        for I in INDICES:

            te = SELF.mesh.trace.elements[I]
            assert te.IS_on_mesh_boundary, f"must be the case because ftype == 'boundary-wise!"

            # only accept 1d arrays
            xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, parse_3_1d_eps=True)

            bn = te.on_mesh_boundary
            assert bn in func, f"trace element #{I} is on <{bn}> which is not covered by boundary-wise func."
            func_i = func[bn]

            vx_i = func_i[0](*xyz_i)
            vy_i = func_i[1](*xyz_i)
            vz_i = func_i[2](*xyz_i)

            if ravel:
                xyz[I] = [_.ravel('F') for _ in xyz_i]
                value[I] = [vx_i.ravel('F'), vy_i.ravel('F'), vz_i.ravel('F')]
            else:
                xyz[I] = xyz_i
                value[I] = [vx_i, vy_i, vz_i, ]


        return xyz, value