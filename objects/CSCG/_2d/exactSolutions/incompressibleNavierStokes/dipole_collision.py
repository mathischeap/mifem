# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 2:37 PM
See [H.J.H. Clercx, C.-H. Bruneau, The normal and oblique collision of a dipole with a no-slip
boundary, Computers & Fluids, 2006]

and

[G.G. de Diego, A. Palha, M. Gerritsma, Inclusion of no-slip boundary conditions in the MEEVC
scheme, JCP, 2019]

"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np

from objects.CSCG._2d.exactSolutions.incompressibleNavierStokes.base import \
    incompressibleNavierStokesBase
from objects.CSCG._2d.fields.vector.main import _2dCSCG_VectorField


class DipoleCollision(incompressibleNavierStokesBase):
    """The 2d dipole collision case. See [H.J.H. Clercx et al., Computers & Fluids, 2006]"""
    def __init__(self, mesh, nu=1/625, x1=0, y1=0.1, x2=0, y2=-0.1, r0=0.1, we1=320, we2=-320):
        """
        (x1, y1) : center of the monopole 1
        (x2, y2) : center of the monopole 2

        These default inputs imply a typical normal collision.

        Parameters
        ----------
        mesh
        nu
        x1
        y1
        x2
        y2
        r0
        we1
        we2
        """
        self._x1_ = x1
        self._y1_ = y1
        self._x2_ = x2
        self._y2_ = y2
        self._r0_ = r0
        self._r0_square_ = r0 ** 2
        self._we1_ = we1
        self._we2_ = we2
        assert we1 + we2  == 0
        self._abs_we_ = abs(we1)
        super(DipoleCollision, self).__init__(mesh, nu)

    @property
    def valid_time(self):
        return "valid_only_at_its_first_instant"

    @property
    def Re(self):
        return 1 / self.nu

    def fx(self, t, x, y):
        return 0 * x

    def fy(self, t, x, y):
        return 0 * y

    @property
    def body_force(self):
        """body force valid all the time."""
        if self._body_force_ is None:
            self._body_force_ =_2dCSCG_VectorField(
                self.mesh, [self.fx, self.fy], valid_time=None, name='body_force')
        return self._body_force_


    def u(self, t, x, y):
        """"""
        x1 = self._x1_
        y1 = self._y1_
        x2 = self._x2_
        y2 = self._y2_

        r1_square = (x - x1) ** 2 + (y - y1) ** 2
        r2_square = (x - x2) ** 2 + (y - y2) ** 2

        term1 = -0.5 * self._abs_we_ * (y - y1) * np.exp(- r1_square / self._r0_square_)
        term2 =  0.5 * self._abs_we_ * (y - y2) * np.exp(- r2_square / self._r0_square_)

        return term1 + term2

    def v(self, t ,x ,y):
        x1 = self._x1_
        y1 = self._y1_
        x2 = self._x2_
        y2 = self._y2_

        r1_square = (x - x1) ** 2 + (y - y1) ** 2
        r2_square = (x - x2) ** 2 + (y - y2) ** 2

        term1 =  0.5 * self._abs_we_ * (x - x1) * np.exp(- r1_square / self._r0_square_)
        term2 = -0.5 * self._abs_we_ * (x - x2) * np.exp(- r2_square / self._r0_square_)

        return term1 + term2

    # def omega(self, t, x, y):
    #     """"""
    #     x1 = self._x1_
    #     y1 = self._y1_
    #     x2 = self._x2_
    #     y2 = self._y2_
    #     we1 = self._we1_
    #     we2 = self._we2_
    #
    #     r1_square = (x - x1) ** 2 + (y - y1) ** 2
    #     r2_square = (x - x2) ** 2 + (y - y2) ** 2
    #
    #
    #     omega_1 = we1 * (1 - r1_square/ self._r0_square_) * np.exp(- r1_square/ self._r0_square_)
    #     omega_2 = we2 * (1 - r2_square/ self._r0_square_) * np.exp(- r2_square/ self._r0_square_)
    #
    #     return omega_2 + omega_1


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/exact_solutions/status/incompressibleNavierStokes/dipole_collision.py
    from objects.CSCG._2d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', bounds=[[-1, 1], [-1, 1]], c=0.0)([2, 2])
    es = ExactSolutionSelector(mesh)("icpsNS:dipole_collision", show_info=True)

    # es.status.vorticity.visualize(time=0)
    a = es.___Pr_compute_Ln_norm_of___('vorticity')
    print(0.5 * a**2 * 2 / 2.282512)