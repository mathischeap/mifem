# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/05 4:14 PM
"""

from components.freeze.base import FrozenOnly


class OnTraceElement_Standard(FrozenOnly):
    """"""
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
        sf = self._sf_

        xyz = dict()
        value = dict()

        func = sf.___DO_evaluate_func_at_time___()

        if i in (None, 'on_mesh_boundaries'):
            RTE = sf.mesh.boundaries.range_of_trace_elements
            INDICES = list()
            for bn in sf.mesh.boundaries.names:
                INDICES.extend(RTE[bn])

        else:
            raise NotImplementedError(f"_3dCSCG_ScalarField of ftype 'standard'"
                                      f"trace-element-reconstruction currently doesn't accept i={i}.")

        for _I in INDICES:
            te = sf.mesh.trace.elements[_I]
            xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, parse_3_1d_eps=True)
            func_i = func[0]

            v_i = func_i(*xyz_i)

            if ravel:
                xyz[_I] = [_.ravel('F') for _ in xyz_i]
                value[_I] = [v_i.ravel('F'), ]
            else:
                xyz[_I] = xyz_i
                value[_I] = [v_i, ]

        return xyz, value
