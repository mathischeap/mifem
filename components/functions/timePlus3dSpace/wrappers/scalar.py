# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/13/2022 2:54 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly

from components.numerical.timePlus3dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txyz_Functions

from importlib import import_module

from components.functions.timePlus3dSpace.wrappers.helpers.scalar_mul import t3d_ScalarMultiply
from components.functions.timePlus3dSpace.wrappers.helpers._3scalars_add import t3d_3ScalarAdd
from components.functions.timePlus3dSpace.wrappers.helpers.scalar_add import t3d_ScalarAdd
from components.functions.timePlus3dSpace.wrappers.helpers.scalar_sub import t3d_ScalarSub
from components.functions.timePlus3dSpace.wrappers.helpers.scalar_neg import t3d_ScalarNeg

from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField


class t3dScalar(FrozenOnly):
    """"""

    def __init__(self, s):
        """"""
        self._s_ = s
        self.__NPD__ = None
        self._freeze_self_()

    def __call__(self, t, x, y, z):
        return self._s_(t, x, y, z)

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
            scalar = _3dCSCG_ScalarField(mesh, self._s_)
            scalar.current_time = t
            return scalar.visualize

        else:
            raise NotImplementedError()

    @property
    def ndim(self):
        return 3

    @property
    def _NPD_(self):
        if self.__NPD__ is None:
            self.__NPD__ = NumericalPartialDerivative_txyz_Functions(self._s_)
        return self.__NPD__

    @property
    def time_derivative(self):
        ps_pt = self._NPD_('t')
        return self.__class__(ps_pt)

    @property
    def gradient(self):
        """"""

        px = self._NPD_('x')
        py = self._NPD_('y')
        pz = self._NPD_('z')

        base_path = '.'.join(str(self.__class__).split(' ')[1][1:].split('.')[:-2]) + '.'
        V_CLASS = getattr(import_module(base_path + "vector"), "t3dVector")

        return V_CLASS(px, py, pz)

    def convection_by(self, u):
        """We compute (u cdot nabla) of self.

        Parameters
        ----------
        u

        Returns
        -------

        """
        assert u.__class__.__name__ == "t3dVector", f"I need a t3dVector."

        px = self._NPD_('x')
        py = self._NPD_('y')
        pz = self._NPD_('z')

        vx, vy, vz = u._v0_, u._v1_, u._v2_

        sx = t3d_ScalarMultiply(vx, px)
        sy = t3d_ScalarMultiply(vy, py)
        sz = t3d_ScalarMultiply(vz, pz)

        return self.__class__(t3d_3ScalarAdd(sx, sy, sz))

    def __add__(self, other):
        """"""
        if other.__class__ is self.__class__:

            s0_add_s1 = t3d_ScalarAdd(self._s_, other._s_)

            return self.__class__(s0_add_s1)

        else:
            raise NotImplementedError()

    def __sub__(self, other):
        """"""
        if other.__class__ is self.__class__:

            s0_sub_s1 = t3d_ScalarSub(self._s_, other._s_)

            return self.__class__(s0_sub_s1)

        else:
            raise NotImplementedError()

    def __neg__(self):
        """"""

        neg = t3d_ScalarNeg(self._s_)

        return self.__class__(neg)

    def __mul__(self, other):
        """"""
        if other.__class__ is self.__class__:
            s0_mul_s1 = t3d_ScalarMultiply(self._s_, other._s_)
            return self.__class__(s0_mul_s1)

        else:
            raise NotImplementedError()
