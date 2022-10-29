# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 3:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.standard._0.base.main import mpRfT2_S0F


class mpRfT2_So0F(mpRfT2_S0F):
    """"""

    def __init__(self, mesh, hybrid=True,
        numbering_parameters='Naive',  name='inner-oriented-0-form'):
        """

        Parameters
        ----------
        mesh
        hybrid :
            {True,}
        numbering_parameters
        name
        """
        super(mpRfT2_So0F, self).__init__(mesh, hybrid, 'outer', numbering_parameters, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_standard_outer_0_form')
        self._freeze_self_()









if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_0/outer/main.py

    from __init__ import rfT2

    mesh = rfT2.rm(500, N_range=(2,3), refinement_intensity=0.2)

    f = mpRfT2_So0F(mesh)

    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar
    import numpy as np
    def p(t, x, y): return np.sin(np.pi * x) * np.cos(np.pi * y) + t
    s = mpRfT2_Scalar(mesh, p)

    f.TW.func = s
    s.current_time = 0

    f.discretization()
    # R = f.reconstruct()
    f.visualization()
    # print(f.error.L())

    df = f.coboundary()

    df.visualization()

    # print(f.num.local_dofs)