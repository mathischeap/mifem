# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/05 1:58 PM

See section 4.1 of [Hu, Lee, Xu, Helicity-conservative finite element discretization for
incompressible MHD systems, JCP, (2021)]

"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.CSCG._3d.exact_solutions.status.incompressible_MHD.base import incompressible_MHD_Base


class AS1(incompressible_MHD_Base):
    """"""

    def __init__(self, es, Re=1e4, Rm=1e4, c=1):
        """"""
        super(AS1, self).__init__(es, Re, Rm, c)
        self._freeze_self_()

    @staticmethod
    def _h(a): return (a ** 2 - a) ** 2
    @staticmethod
    def _h_p(a): return 2 * (a ** 2 - a)* (2*a - 1)
    @staticmethod
    def _g1(t): return 4 - 2 * t
    @staticmethod
    def _g2(t): return 1 + t
    @staticmethod
    def _g3(t): return 1 - t

    def u(self, t, x, y, z):
        return - self._g1(t) * self._h_p(x) * self._h(y) * self._h(z)
    def v(self, t, x, y, z):
        return - self._g2(t) * self._h(x) * self._h_p(y) * self._h(z)
    def w(self, t, x, y, z):
        return - self._g3(t) * self._h(x) * self._h(y) * self._h_p(z)

    # def u_y(self, t, x, y, z):
    #     return - self._g1(t) * self._h_p(x) * self._h_p(y) * self._h(z)
    # def u_z(self, t, x, y, z):
    #     return - self._g1(t) * self._h_p(x) * self._h(y) * self._h_p(z)
    #
    # def v_x(self, t, x, y, z):
    #     return - self._g2(t) * self._h_p(x) * self._h_p(y) * self._h(z)
    # def v_z(self, t, x, y, z):
    #     return - self._g2(t) * self._h(x) * self._h_p(y) * self._h_p(z)
    #
    # def w_x(self, t, x, y, z):
    #     return - self._g3(t) * self._h_p(x) * self._h(y) * self._h_p(z)
    # def w_y(self, t, x, y, z):
    #     return - self._g3(t) * self._h(x) * self._h_p(y) * self._h_p(z)

    def Bx(self, t, x, y, z):
        return self.w_y(t, x, y, z) - self.v_z(t, x, y, z)
    def By(self, t, x, y, z):
        return self.u_z(t, x, y, z) - self.w_x(t, x, y, z)
    def Bz(self, t, x, y, z):
        return self.v_x(t, x, y, z) - self.u_y(t, x, y, z)

    def s(self, t, x, y, z):
        """Non-zero mass source term."""
        return self.u_x(t, x, y, z) + self.v_y(t, x, y, z) + self.w_z(t, x, y, z)

    def p(self, t, x, y, z):
        """"""
        return self._h(x) * self._h(y) * self._h(z)





if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_3d/exact_solutions/status/incompressible_MHD/analytic_solution_1.py

    from objects.CSCG._3d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', c=0.0)([5, 5, 5])
    es = ExactSolutionSelector(mesh)("MHD:as1", show_info=True)

    r = es.status.volume_current_density
    r.current_time = 1
    r.visualize()