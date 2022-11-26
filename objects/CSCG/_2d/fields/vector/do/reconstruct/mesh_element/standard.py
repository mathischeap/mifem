# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly





class OnMeshElement_for_Standard(FrozenOnly):
    """"""
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self,xi, eta, ravel, i):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param i:
        :return:
        """
        SELF = self._vf_

        xyz = dict()
        value = dict()

        if isinstance(i, int):
            INDICES = [i, ]
        elif isinstance(i, (tuple, list)):
            INDICES = i
        elif i is None:
            INDICES = SELF.mesh.elements.indices
        else:
            raise NotImplementedError(f"reconstruction do not accept i={i}.")

        func = SELF.___DO_evaluate_func_at_time___()
        for i in INDICES:
            element = SELF.mesh.elements[i]
            xyz_i = element.coordinate_transformation.mapping(xi, eta)
            vx_i = func[0](*xyz_i)
            vy_i = func[1](*xyz_i)

            if ravel:
                xyz[i] = [I.ravel('F') for I in xyz_i]
                value[i] = [vx_i.ravel('F'), vy_i.ravel('F')]
            else:
                xyz[i] = xyz_i
                value[i] = [vx_i, vy_i]

        return xyz, value