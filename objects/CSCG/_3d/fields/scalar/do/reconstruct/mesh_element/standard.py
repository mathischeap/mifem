# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly


class OnMeshElement_Standard(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel, i):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param i:
            In which mesh-element(s) we reconstruct the scalar.

            if `i` is:
                1) int: reconstruct in the particular mesh element #i.
                2) None: in all mesh-elements.

        :return:
        """
        SELF = self._sf_

        xyz = dict()
        value = dict()

        if isinstance(i, int):
            INDICES = [i,]
        elif i is None:
            INDICES = SELF.mesh.elements.indices
        else:
            raise NotImplementedError( f"We can not reconstruct in elements={i}.")

        func = SELF.___DO_evaluate_func_at_time___()

        for i in INDICES:
            element = SELF.mesh.elements[i]
            xyz_i = element.coordinate_transformation.mapping(xi, eta, sigma)
            v_i = func[0](*xyz_i)

            if ravel:
                xyz[i] = [I.ravel('F') for I in xyz_i]
                value[i] = [v_i.ravel('F'),]
            else:
                xyz[i] = xyz_i
                value[i] = [v_i,]

        return xyz, value