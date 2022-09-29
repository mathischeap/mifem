# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 2:14 PM
"""
from screws.freeze.base import FrozenOnly


class OnMeshElement_for_Standard(FrozenOnly):
    """"""
    def __init__(self, scalar):
        self._scalar_ = scalar
        self._freeze_self_()

    def __call__(self, xi, eta, ravel, i):
        """

        :param xi: We will map (xi, eta) to physical domain without any processing like `meshgrid`.
        :param eta: We will map (xi, eta) to physical domain without any processing like `meshgrid`.
        :param ravel:
        :param i:
            1) self.ftype == 'standard':
                Do the reconstruction in mesh element #i. When it is None, it means all local mesh
                elements.
        :return:
        """
        SELF = self._scalar_

        xyz = dict()
        value = dict()

        # parse elements -----------------------------------------------------------------------
        if isinstance(i, int):
            INDICES = [i, ]
        elif isinstance(i, (tuple, list)):
            INDICES = i
        elif i is None:
            INDICES = SELF.mesh.elements.indices
        else:
            raise NotImplementedError(f"reconstruction do not accept i={i}.")

        #------- do the reconstruction ------------------------------------------------------
        func = SELF.___DO_evaluate_func_at_time___()[0]
        for i in INDICES:
            element = SELF.mesh.elements[i]
            xyz_i = element.coordinate_transformation.mapping(xi, eta)
            v_i = func(*xyz_i)

            if ravel:
                xyz[i] = [I.ravel('F') for I in xyz_i]
                value[i] = [v_i.ravel('F'), ]
            else:
                xyz[i] = xyz_i
                value[i] = [v_i, ]

        return xyz, value