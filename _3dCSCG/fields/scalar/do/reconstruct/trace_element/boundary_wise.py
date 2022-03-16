from screws.freeze.inheriting.frozen_only import FrozenOnly


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
            1) self.ftype == 'standard':
                Do the reconstruction in mesh element #i. When it is None, it means all local mesh
                elements.
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

        for I in INDICES:
            te = SELF.mesh.trace.elements[I]
            assert te.IS_on_mesh_boundary, f"must be the case!"
            # when we do trace-element-wise reconstruction, we only accept 1d xi, eta, sigma.
            xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, parse_3_1d_eps=True)

            bn = te.on_mesh_boundary
            assert bn in func, f"trace element #{I} is on <{bn}> which is not covered by boundary-wise func."
            func_i = func[bn][0]

            v_i = func_i(*xyz_i)

            if ravel:
                xyz[I] = [_.ravel('F') for _ in xyz_i]
                value[I] = [v_i.ravel('F'), ]
            else:
                xyz[I] = xyz_i
                value[I] = [v_i, ]

        return xyz, value