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
            1) self.ftype == 'standard':
                Do the reconstruction in mesh element #i. When it is None, it means all local mesh
                elements.
        :return:
        """
        SELF = self._sf_

        xyz = dict()
        value = dict()

        assert isinstance(i, int) or i is None, f"We currently only accept int or None for i"
        INDICES = SELF.mesh.elements.indices if i is None else [i, ]

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