# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/26 2:19 PM
"""

import numpy as np

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.trace.elements.element.coordinate_transformation.constant import _3dTraceElement_CT_Constant


class _3dCSCG_Trace_Element_CoordinateTransformation(FrozenOnly):
    def __init__(self, te):
        self._te_ = te
        self._constant_ = _3dTraceElement_CT_Constant(te)
        self._freeze_self_()

    @property
    def constant(self):
        return self._constant_

    def mapping(self, *evaluation_points, from_element=None, side=None, parse_3_1d_eps=False):
        """
        The local mapping.

        :param evaluation_points : A tuple or list of shape (ndim-1, ...).
        :param from_element: (default: ``None``) We try to compute the
            mapping from a given mesh element. When it is None,
            we compute it from its CHARACTERISTIC_element. Notice that
            when a trace element in on periodic domain, we have to
            provide a ``from_element``. Because this is the case,
            we will get different mapping from different element.
        :param side: when one mesh element is periodic to itself, we may
            need to provide side as well.
        :param parse_3_1d_eps: If `parse_ep` is True, then we have *ep is xi, eta, sigma, and they are all 1d,
            between [-1,1], we will pick up two from them according the side and do the mesh grid.
        """
        if self._te_.IS.on_periodic_boundary:
            assert from_element is not None, \
                "to compute the mapping of a trace element on periodic " \
                "boundary, we have to provide from which element you " \
                "want to compute it since it clearly will gives " \
                "different results."
            if from_element == 'any':
                # the different results do not matter; for example, when
                # we want to get value from a periodic function, the
                # location for evaluating the function also does not
                # matter.
                from_element = self._te_.CHARACTERISTIC_element
            else:
                pass

        if from_element is None:
            i = self._te_.CHARACTERISTIC_element
        elif from_element == 'any':
            i = self._te_.CHARACTERISTIC_element
        else:
            i = from_element

        assert self._te_.i in self._te_._elements_.map[i], \
            f"trace element {self._te_.i} is not on mesh element {i}."

        if self._te_._elements_.map[i].count(self._te_.i) == 1:
            side_index = self._te_._elements_.map[i].index(self._te_.i)
            element_side = 'NSWEBF'[side_index]
            if from_element is None:
                assert element_side == self._te_.CHARACTERISTIC_side

            if side is not None:
                assert element_side == side, f"cannot compute on provided side {side}"

        elif self._te_._elements_.map[i].count(self._te_.i) == 2:
            assert side is not None, f"trace element #{self._te_.i} " \
                                     f"is on two sides of element #{i} " \
                                     f"(periodic), provide side as well."
            element_side = side
        else:
            raise Exception()

        assert self._te_.i == self._te_._elements_.map[i]['NSWEBF'.index(element_side)], \
            f"trace element #{self._te_.i} is not at {element_side} of mesh element #{i}."

        ep = self._te_._elements_.___generate_full_ep___(evaluation_points,
                                                         element_side,
                                                         parse_3_1d_eps=parse_3_1d_eps)
        x, y, z = self._te_._mesh_.elements[i].coordinate_transformation.mapping(*ep)
        return x, y, z

    def Jacobian_matrix(self, *evaluation_points, parse_3_1d_eps=False):
        """
        The local Jacobian matrix.

        :param evaluation_points: A tuple or list of shape (ndim-1, ...).
        :param parse_3_1d_eps:
        """
        i = self._te_.CHARACTERISTIC_element
        element_side = self._te_.CHARACTERISTIC_side
        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_side, parse_3_1d_eps=parse_3_1d_eps)
        J = self._te_._mesh_.elements[i].coordinate_transformation.Jacobian_matrix(*ep)
        if element_side in 'NS':
            return ((J[0][1], J[0][2]),
                    (J[1][1], J[1][2]),
                    (J[2][1], J[2][2]))
        elif element_side in 'WE': # this is very important, do not used [0], [2] for the indices. But when we feed x, y, z, always use x, y, z.
            return ((J[0][2], J[0][0]),
                    (J[1][2], J[1][0]),
                    (J[2][2], J[2][0]))
        else:
            return ((J[0][0], J[0][1]),
                    (J[1][0], J[1][1]),
                    (J[2][0], J[2][1]))

    def inverse_Jacobian_matrix(self, *evaluation_points):
        """
        The local inverse_Jacobian matrix.

        :param evaluation_points : A tuple or list of shape (ndim-1, ...).
        """
        i = self._te_.CHARACTERISTIC_element
        element_side = self._te_.CHARACTERISTIC_side
        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_side)
        iJ = self._te_._mesh_.elements[i].coordinate_transformation.inverse_Jacobian_matrix(*ep)
        if element_side in 'NS':
            return ((iJ[1][0], iJ[1][1], iJ[1][2]),
                    (iJ[2][0], iJ[2][1], iJ[2][2]))
        elif element_side in 'WE': # this is very important, do not used [0], [2] for the indices.. But when we feed x, y, z, always use x, y, z.
            return ((iJ[2][0], iJ[2][1], iJ[2][2]),
                    (iJ[0][0], iJ[0][1], iJ[0][2]))
        else:
            return ((iJ[0][0], iJ[0][1], iJ[0][2]),
                    (iJ[1][0], iJ[1][1], iJ[1][2]))

    def metric_matrix(self, *evaluation_points):
        """Compute the metric matrix G whose entries are g_{i,j}."""
        J = self.Jacobian_matrix(*evaluation_points)
        Gk = [[None for _ in range(2)] for __ in range(2)]
        for i in range(2):
            for j in range(i, 2):
                Gk[i][j] = J[0][i] * J[0][j]
                for l in range(1, 3):
                    Gk[i][j] += J[l][i] * J[l][j]
                if i != j:
                    Gk[j][i] = Gk[i][j]
        return Gk


    def inverse_metric_matrix(self, *evaluation_points):
        """Compute the inverse metric matrix G whose entries are g^{i,j}."""
        iJ = self.inverse_Jacobian_matrix(*evaluation_points)
        Gk = [[None for _ in range(2)] for __ in range(2)]
        for i in range(2):
            for j in range(i, 2):
                Gk[i][j] = iJ[i][0] * iJ[j][0] + iJ[i][1] * iJ[j][1] + iJ[i][2] * iJ[j][2]
                if i != j:
                    Gk[j][i] = Gk[i][j]

        return Gk



    def metric(self, *evaluation_points):
        """return metric g."""
        G = self.metric_matrix(*evaluation_points)
        # noinspection PyUnresolvedReferences
        return G[0][0]*G[1][1] - G[0][1]*G[1][0]

    def ___PRIVATE_outward_unit_normal_vector___(self, *evaluation_points, from_element=None, side=None):
        """For a single trace element (a surface in 3D space), it is
        hard to say which direction is the positive direction (the
        outward direction). So we put onto the surface of the mesh
        element. So, this outward unit norm vector of this trace element
        is the outward unit norm vector of the characteristic mesh
        element on the face the trace element is located.

        In fact, for a trace element, it seems not good to consider it
        has a normal vector. We should consider so for a mesh element.
        For example, when we have a 2-form vector(u), the trace of it,
        vector(u) dot vector(n), gives a scalar field anyway, the
        outward normal vector (n) is used at the mesh element level, as
        well as the trace matrix N.

        vec(a) = (a1, a2, a3)
        vec(b) = (b1, b2, b3)

        vec(a) x vec(b) = (a2 b3 - a3 b2, a3 b1 - a1 b3, a1 b2-a2 b1)

        :param evaluation_points: A tuple or list of shape (ndim-1, ...).
        :param from_element: We will compute the normal vector by
            considering it is on this mesh element. Please give this
            parameter in case of random error!
        :param side: when one mesh element is periodic to itself, we may
            need to provide side as well.
        """
        J = self.Jacobian_matrix(*evaluation_points)
        a = (J[0][0], J[1][0], J[2][0])
        b = (J[0][1], J[1][1], J[2][1])
        acb0 = a[1] * b[2] - a[2] * b[1]
        acb1 = a[2] * b[0] - a[0] * b[2]
        acb2 = a[0] * b[1] - a[1] * b[0]
        norm = np.sqrt(acb0**2 + acb1**2 + acb2**2)

        if from_element is None:
            i = self._te_.CHARACTERISTIC_element
        else:
            i = from_element
        assert self._te_.i in self._te_._elements_.map[i], \
            f"trace element {self._te_.i} is not on mesh element {i}."

        if self._te_._elements_.map[i].count(self._te_.i) == 1:
            side_index = self._te_._elements_.map[i].index(self._te_.i)
            element_side = 'NSWEBF'[side_index]
            if from_element is None:
                assert element_side == self._te_.CHARACTERISTIC_side

            if side is not None:
                assert element_side == side, f"cannot compute on provided side {side}"

        elif self._te_._elements_.map[i].count(self._te_.i) == 2:
            assert side is not None, f"trace element #{self._te_.i} " \
                                     f"is on two sides of element #{i} " \
                                     f"(periodic), provide side as well."
            element_side = side
        else:
            raise Exception()

        if side is not None: element_side = side #

        uv = np.array([acb0 / norm, acb1 / norm, acb2 / norm])
        if element_side in 'NWB':
            uv = - uv
        else:
            pass

        return uv

    def unit_normal_vector(self, *evaluation_points, parse_3_1d_eps=False):
        """The unit_normal_vector (vector n) according to the right-hand-rule (without considering it is
        which side of the mesh element; it can point inner direction of the mesh element).

        vec(a) = (a1, a2, a3)
        vec(b) = (b1, b2, b3)

        vec(a) x vec(b) = (a2 b3 - a3 b2, a3 b1 - a1 b3, a1 b2-a2 b1)

        :param evaluation_points: A tuple or list of shape (ndim-1, ...).
        :param parse_3_1d_eps:
        """
        J = self.Jacobian_matrix(*evaluation_points, parse_3_1d_eps=parse_3_1d_eps)
        a = (J[0][0], J[1][0], J[2][0])
        b = (J[0][1], J[1][1], J[2][1])
        acb0 = a[1] * b[2] - a[2] * b[1]
        acb1 = a[2] * b[0] - a[0] * b[2]
        acb2 = a[0] * b[1] - a[1] * b[0]
        norm = np.sqrt(acb0**2 + acb1**2 + acb2**2)

        return np.array([acb0 / norm, acb1 / norm, acb2 / norm])
