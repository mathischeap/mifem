# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/3/2022 6:29 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly

from components.numerical.timePlus3dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txyz_Functions

from importlib import import_module
from components.functions.timePlus3dSpace.wrappers.helpers._3scalars_add import t3d_3ScalarAdd
from components.functions.timePlus3dSpace.wrappers.helpers.scalar_sub import t3d_ScalarSub
from components.functions.timePlus3dSpace.wrappers.helpers.scalar_mul import t3d_ScalarMultiply
from components.functions.timePlus3dSpace.wrappers.helpers.scalar_add import t3d_ScalarAdd
from components.functions.timePlus3dSpace.wrappers.helpers.scalar_neg import t3d_ScalarNeg


from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField


class t3dVector(FrozenOnly):
    """ Wrap three functions into a vector class.
    """

    def __init__(self, v0, v1, v2):
        """Initialize a vector with 3 functions which take (t, x, y, z) as inputs.

        Parameters
        ----------
        v0
        v1
        v2
        """
        self._v0_ = v0
        self._v1_ = v1
        self._v2_ = v2
        self._vs_ = [v0, v1, v2]
        self.__NPD0__ = None
        self.__NPD1__ = None
        self.__NPD2__ = None
        self._freeze_self_()

    def __call__(self, t, x, y, z):
        """Evaluate the vector at (t, x, y, z)"""
        return self._v0_(t, x, y, z), self._v1_(t, x, y, z), self._v2_(t, x, y, z)

    def __getitem__(self, item):
        return self._vs_[item]

    def visualize(self, mesh, t):
        """Return a visualize class for a mesh at t=`t`.

        Parameters
        ----------
        mesh
        t

        Returns
        -------

        """
        if mesh.__class__.__name__ == "_3dCSCG_Mesh":
            vector = _3dCSCG_VectorField(mesh, self._vs_)
            vector.current_time = t
            return vector.visualize

        else:
            raise NotImplementedError()

    @property
    def ndim(self):
        return 3

    @property
    def _NPD0_(self):
        if self.__NPD0__ is None:
            self.__NPD0__ = NumericalPartialDerivative_txyz_Functions(self._v0_)
        return self.__NPD0__

    @property
    def _NPD1_(self):
        if self.__NPD1__ is None:
            self.__NPD1__ = NumericalPartialDerivative_txyz_Functions(self._v1_)
        return self.__NPD1__

    @property
    def _NPD2_(self):
        if self.__NPD2__ is None:
            self.__NPD2__ = NumericalPartialDerivative_txyz_Functions(self._v2_)
        return self.__NPD2__

    @property
    def time_derivative(self):
        pv0_pt = self._NPD0_('t')
        pv1_pt = self._NPD1_('t')
        pv2_pt = self._NPD2_('t')
        return self.__class__(pv0_pt, pv1_pt, pv2_pt)

    @property
    def divergence(self):
        pv0_px = self._NPD0_('x')
        pv1_py = self._NPD1_('y')
        pv2_pz = self._NPD2_('z')

        scalar_function = t3d_3ScalarAdd(pv0_px, pv1_py, pv2_pz)

        base_path = '.'.join(str(self).split(' ')[0][1:].split('.')[:-2]) + '.'
        S_CLASS = getattr(import_module(base_path + "scalar"), "t3dScalar")

        return S_CLASS(scalar_function)

    @property
    def curl(self):
        """The curl of a 3d vector. Lets say self is (u, v, w):

        (pw/py - pv/pz, pu/pz - pw/px, pv/px - pu/py)

        Returns
        -------

        """
        pw_py = self._NPD2_('y')
        pw_px = self._NPD2_('x')
        pu_pz = self._NPD0_('z')
        pu_py = self._NPD0_('y')
        pv_px = self._NPD1_('x')
        pv_pz = self._NPD1_('z')

        v0 = t3d_ScalarSub(pw_py, pv_pz)
        v1 = t3d_ScalarSub(pu_pz, pw_px)
        v2 = t3d_ScalarSub(pv_px, pu_py)

        return self.__class__(v0, v1, v2)

    def convection_by(self, u):
        """We compute (u cdot nabla) of self where u is another t3d vector.

        Parameters
        ----------
        u

        Returns
        -------

        """
        assert u.__class__.__name__ == "t3dVector", f"I need a t3dVector."

        v0px = self._NPD0_('x')
        v0py = self._NPD0_('y')
        v0pz = self._NPD0_('z')
        v1px = self._NPD1_('x')
        v1py = self._NPD1_('y')
        v1pz = self._NPD1_('z')
        v2px = self._NPD2_('x')
        v2py = self._NPD2_('y')
        v2pz = self._NPD2_('z')

        vx, vy, vz = u._v0_, u._v1_, u._v2_

        v0x = t3d_ScalarMultiply(vx, v0px)
        v0y = t3d_ScalarMultiply(vy, v0py)
        v0z = t3d_ScalarMultiply(vz, v0pz)
        v1x = t3d_ScalarMultiply(vx, v1px)
        v1y = t3d_ScalarMultiply(vy, v1py)
        v1z = t3d_ScalarMultiply(vz, v1pz)
        v2x = t3d_ScalarMultiply(vx, v2px)
        v2y = t3d_ScalarMultiply(vy, v2py)
        v2z = t3d_ScalarMultiply(vz, v2pz)

        return self.__class__(t3d_3ScalarAdd(v0x, v0y, v0z),
                              t3d_3ScalarAdd(v1x, v1y, v1z),
                              t3d_3ScalarAdd(v2x, v2y, v2z))

    def __add__(self, other):
        """

        Parameters
        ----------
        other

        Returns
        -------

        """
        if other.__class__ is self.__class__:

            v00, v01, v02 = self._v0_, self._v1_, self._v2_
            v10, v11, v12 = other._v0_, other._v1_, other._v2_

            V0 = t3d_ScalarAdd(v00, v10)
            V1 = t3d_ScalarAdd(v01, v11)
            V2 = t3d_ScalarAdd(v02, v12)

            return self.__class__(V0, V1, V2)

        else:
            raise NotImplementedError()

    def __sub__(self, other):
        """

        Parameters
        ----------
        other

        Returns
        -------

        """
        if other.__class__ is self.__class__:

            v00, v01, v02 = self._v0_, self._v1_, self._v2_
            v10, v11, v12 = other._v0_, other._v1_, other._v2_

            V0 = t3d_ScalarSub(v00, v10)
            V1 = t3d_ScalarSub(v01, v11)
            V2 = t3d_ScalarSub(v02, v12)

            return self.__class__(V0, V1, V2)

        else:
            raise NotImplementedError()

    def __neg__(self):
        v0, v1, v2 = self._v0_, self._v1_, self._v2_

        neg_v0 = t3d_ScalarNeg(v0)
        neg_v1 = t3d_ScalarNeg(v1)
        neg_v2 = t3d_ScalarNeg(v2)

        return self.__class__(neg_v0, neg_v1, neg_v2)

    def dot(self, other):
        """`self` dot product with `other`. So lets say self = (a, b), other = (u, v),
        self.dot(other) gives a scalar, au + bv.

        Parameters
        ----------
        other

        Returns
        -------

        """
        if other.__class__ is self.__class__:

            v00, v01, v02 = self._v0_, self._v1_, self._v2_
            v10, v11, v12 = other._v0_, other._v1_, other._v2_

            V0 = t3d_ScalarMultiply(v00, v10)
            V1 = t3d_ScalarMultiply(v01, v11)
            V2 = t3d_ScalarMultiply(v02, v12)

            V0V1V2 = t3d_3ScalarAdd(V0, V1, V2)
            base_path = '.'.join(str(self).split(' ')[0][1:].split('.')[:-2]) + '.'
            S_CLASS = getattr(import_module(base_path + "scalar"), "t3dScalar")
            return S_CLASS(V0V1V2)

        else:
            raise NotImplementedError()


    def cross_product(self, other):
        """"""
        if other.__class__ is self.__class__:

            a, b, c = self._v0_, self._v1_, self._v2_
            u, v, w = other._v0_, other._v1_, other._v2_

            bw = t3d_ScalarMultiply(b, w)
            cu = t3d_ScalarMultiply(c, u)
            av = t3d_ScalarMultiply(a, v)
            cv = t3d_ScalarMultiply(c, v)
            aw = t3d_ScalarMultiply(a, w)
            bu = t3d_ScalarMultiply(b, u)

            V0 = self.__class__(bw, cu, av)
            V1 = self.__class__(cv, aw, bu)

            return V0 - V1

        else:
            raise NotImplementedError()
