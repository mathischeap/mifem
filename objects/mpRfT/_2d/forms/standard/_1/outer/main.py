# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 12:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.standard._1.base.main import mpRfT2_S1F
from objects.mpRfT._2d.forms.standard._1.outer.discretize.main import mpRfT2_So1F_Discretize
from objects.mpRfT._2d.forms.standard._1.outer.reconstruct import mpRfT2_So1F_Reconstruct
from objects.mpRfT._2d.forms.standard._1.outer.migrate import mpRfT2_So1F_Migrate
from objects.mpRfT._2d.forms.standard._1.outer.boundary_integrate import mpRfT2_So1F_BI




class mpRfT2_So1F(mpRfT2_S1F):
    """"""

    def __init__(self, mesh, hybrid=True,
        numbering_parameters='Naive',  name='outer-oriented-1-form'):
        """

        Parameters
        ----------
        mesh
        hybrid :
            {True,}
        numbering_parameters
        name
        """
        super(mpRfT2_So1F, self).__init__(mesh, hybrid, 'outer', numbering_parameters, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_standard_outer_1_form')

        self._discretize_ = mpRfT2_So1F_Discretize(self)
        self._reconstruct_ = mpRfT2_So1F_Reconstruct(self)
        self._migrate_ = mpRfT2_So1F_Migrate(self)
        self._BI_ = mpRfT2_So1F_BI(self)
        self._freeze_self_()


    @property
    def boundary_integrate(self):
        """"""
        return self._BI_








if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/outer/main.py
    from __init__ import rfT2

    fc = rfT2.rf(100, N_range=(3,3))

    f = fc('1-f-o')
    t = fc('nst')

    import numpy as np
    def p(t, x, y): return np.sin(np.pi*x) * np.cos(np.pi*y) + t
    def q(t, x, y): return np.cos(np.pi*x) * np.sin(np.pi*y) + t

    def h(t, x, y): return np.sin(np.pi*x) * np.sin(np.pi*y) + t

    s = fc('scalar', h)
    v = fc('vector', (p, q))

    f.TW.func = v
    v.current_time = 0
    # f.discretize()

    t.TW.func = s
    s.current_time = 0
    # t.discretize()

    # f.boundary_integrate(t)

    print(f.num.local_dofs)