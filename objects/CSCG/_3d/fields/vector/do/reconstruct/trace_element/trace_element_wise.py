# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class OnTraceElement_TraceElementWise(FrozenOnly):
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

        if i is None:  # we reconstruct on all valid local trace elements
            INDICES = list()
            # noinspection PyUnresolvedReferences
            INDICES.extend(func.keys())
        elif i == 'on_mesh_boundaries':
            # we only reconstruct on all the valid local trace elements which are also on mesh boundaries.
            CMB = SELF.covered_mesh_boundaries  # will contain all mesh boundary names.
            RTE = SELF.mesh.boundaries.range_of_trace_elements
            boundary_trace_elements = list()  # local trace elements on all mesh boundaries
            for mb in CMB:
                boundary_trace_elements.extend(RTE[mb])
            ___ = list()
            # noinspection PyUnresolvedReferences
            ___.extend(func.keys())
            INDICES = list()
            for _I in ___:
                if _I in boundary_trace_elements:
                    INDICES.append(_I)

        else:
            raise NotImplementedError(f"_3dCSCG_VectorField of 'trace-element-wise' ftype "
                                      f"trace-element-reconstruction currently don't accept i={i}."
                                      f"i must be one of (None, 'on_mesh_boundaries').")

        for _I in INDICES:  # go through all valid local trace elements

            xyz_i, v_i = func[_I](xi, eta, sigma)

            if ravel:
                xyz[_I] = [_.ravel('F') for _ in xyz_i]
                value[_I] = [_.ravel('F') for _ in v_i]
            else:
                xyz[_I] = xyz_i
                value[_I] = v_i

        return xyz, value
