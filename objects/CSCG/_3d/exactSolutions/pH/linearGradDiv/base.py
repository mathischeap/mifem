# -*- coding: utf-8 -*-
"""
S is a scalar field, V is a vector field. And
dV/dt = - grad S
dS/dt = - div V

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/23/2022 5:35 PM
"""
import sys
import numpy as np

if './' not in sys.path: sys.path.append('./')
from objects.CSCG._3d.exactSolutions.base import Base

from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField

from components.numerical.timePlus3dSpace.partial_derivative_as_functions import NumericalPartialDerivative_txyz_Functions


class pH_LinearGradDiv_Base(Base):
    """"""

    def __init__(self, mesh):
        """"""
        super(pH_LinearGradDiv_Base, self).__init__(mesh)
        self._SF_ = None
        self._VF_ = None
        self._freeze_self_()

    def p(self, t, x, y, z):
        raise NotImplementedError()

    @property
    def sf(self):
        """"""
        if self._SF_ is None:
            self._SF_ = _3dCSCG_ScalarField(
                self.mesh, self.p,
                ftype='standard',
                valid_time=None,
                name='scalar-field')
        return self._SF_


    def u(self, t, x, y, z):
        raise NotImplementedError()

    def v(self, t, x, y, z):
        raise NotImplementedError()

    def w(self, t, x, y, z):
        raise NotImplementedError()

    @property
    def vf(self):
        """"""
        if self._VF_ is None:
            self._VF_ = _3dCSCG_VectorField(
                self.mesh,
                [
                    self.u,
                    self.v,
                    self.w
                ],
                ftype='standard',
                valid_time=None,
                name='vector-field'
            )
        return self._VF_



    def ___PreFrozenChecker___(self):
        """We use this general method to do the check, in particular exact solution, we can define particular
        check method by override this method.
        """
        TS = self.___PRIVATE_generate_random_valid_time_instances___()
        x, y, z = self._mesh_.do.generate_random_coordinates()

        if len(x) == 0: return

        Np = NumericalPartialDerivative_txyz_Functions(self.p)

        px = Np("x")
        py = Np("y")
        pz = Np("z")
        pt = Np("t")

        Nu = NumericalPartialDerivative_txyz_Functions(self.u)
        Nv = NumericalPartialDerivative_txyz_Functions(self.v)
        Nw = NumericalPartialDerivative_txyz_Functions(self.w)

        ut = Nu('t')
        vt = Nv('t')
        wt = Nw('t')

        ux = Nu('x')
        vy = Nv('y')
        wz = Nw('z')

        for t in TS:

            t = float(t)

            np.testing.assert_array_almost_equal(ut(t, x, y, z), - px(t, x, y, z))
            np.testing.assert_array_almost_equal(vt(t, x, y, z), - py(t, x, y, z))
            np.testing.assert_array_almost_equal(wt(t, x, y, z), - pz(t, x, y, z))

            np.testing.assert_array_almost_equal(
                - pt(t, x, y, z),
                ux(t, x, y, z) + vy(t, x, y, z) + wz(t, x, y, z)
            )




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
