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

from components.numerical.timePlus2dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txy_Functions

from components.functions.timePlus2dSpace.wrappers.helpers.scalar_add import t2d_ScalarAdd
from components.functions.timePlus2dSpace.wrappers.helpers.scalar_sub import t2d_ScalarSub
from components.functions.timePlus2dSpace.wrappers.helpers.scalar_neg import t2d_ScalarNeg
from components.functions.timePlus2dSpace.wrappers.helpers.scalar_mul import t2d_ScalarMultiply

from importlib import import_module

from objects.CSCG._2d.fields.scalar.main import _2dCSCG_ScalarField
from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar


class t2dScalar(FrozenOnly):
    """"""

    def __init__(self, s):
        """"""
        self._s_ = s
        self.__NPD__ = None
        self._freeze_self_()

    def __call__(self, t, x, y):
        return self._s_(t, x, y)

    def visualize(self, mesh, t):
        """Return a visualize class for a mesh at t=`t`.

        Parameters
        ----------
        mesh
        t

        Returns
        -------

        """
        if mesh.__class__.__name__ == "_2dCSCG_Mesh":
            scalar = _2dCSCG_ScalarField(mesh, self)
            scalar.current_time = t
            return scalar.visualize

        elif mesh.__class__.__name__ == "miUsGrid_TriangularMesh":
            scalar = miUsGrid_Triangular_Scalar(mesh, self)
            scalar.current_time = t
            return scalar.export
        else:
            raise NotImplementedError()


    @property
    def ndim(self):
        """"""
        return 2

    @property
    def _NPD_(self):
        """"""
        if self.__NPD__ is None:
            self.__NPD__ = NumericalPartialDerivative_txy_Functions(self._s_)
        return self.__NPD__

    @property
    def time_derivative(self):
        """"""
        ps_pt = self._NPD_('t')
        return self.__class__(ps_pt)


    @property
    def gradient(self):
        """"""

        px = self._NPD_('x')
        py = self._NPD_('y')

        base_path = '.'.join(str(self.__class__).split(' ')[1][1:].split('.')[:-2]) + '.'
        V_CLASS = getattr(import_module(base_path + "vector"), "t2dVector")

        return V_CLASS(px, py)
    
    
    @property
    def curl(self):
        """"""
        px = self._NPD_('x')
        py = self._NPD_('y')

        neg_px = - self.__class__(px)

        base_path = '.'.join(str(self.__class__).split(' ')[1][1:].split('.')[:-2]) + '.'
        V_CLASS = getattr(import_module(base_path + "vector"), "t2dVector")

        return V_CLASS(py, neg_px)

    def convection_by(self, u):
        """We compute (u cdot nabla) of self.

        Parameters
        ----------
        u

        Returns
        -------

        """
        assert u.__class__.__name__ == "t2dVector", f"I need a t2dVector."

        px = self._NPD_('x')
        py = self._NPD_('y')

        vx, vy = u._v0_, u._v1_

        sx = t2d_ScalarMultiply(vx, px)
        sy = t2d_ScalarMultiply(vy, py)

        return self.__class__(t2d_ScalarAdd(sx, sy))


    def __add__(self, other):
        """"""
        if other.__class__ is self.__class__:

            s0_add_s1 = t2d_ScalarAdd(self._s_, other._s_)

            return self.__class__(s0_add_s1)

        else:
            raise NotImplementedError()

    def __sub__(self, other):
        """"""
        if other.__class__ is self.__class__:

            s0_sub_s1 = t2d_ScalarSub(self._s_, other._s_)

            return self.__class__(s0_sub_s1)

        else:
            raise NotImplementedError()


    def __neg__(self):
        """"""

        neg = t2d_ScalarNeg(self._s_)

        return self.__class__(neg)


    def __mul__(self, other):
        """"""
        if other.__class__ is self.__class__:
            s0_mul_s1 = t2d_ScalarMultiply(self._s_, other._s_)
            return self.__class__(s0_mul_s1)

        else:
            raise NotImplementedError()


if __name__ == '__main__':
    # mpiexec -n 4 python components/functions/timePlus2dSpace/wrappers/scalar.py
    def f0(t, x, y): return x + y + t
    def f1(t, x, y): return x * y + t

    s0 = t2dScalar(f0)
    s1 = t2dScalar(f1)

    from objects.CSCG._2d.__init__ import mesh

    mesh = mesh('crazy', c=0.1)([5, 5])

    vis0 = s0.visualize(mesh, t=0)
    vis0()
