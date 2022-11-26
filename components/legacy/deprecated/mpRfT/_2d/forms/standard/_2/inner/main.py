# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 12:35 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.forms.standard._2.base.main import mpRfT2_S2F


class mpRfT2_Si2F(mpRfT2_S2F):
    """"""

    def __init__(self, mesh, hybrid=True,
        numbering_parameters='Naive',  name='inner-oriented-2-form'):
        """

        Parameters
        ----------
        mesh
        hybrid :
            {True,}
        numbering_parameters
        name
        """
        super(mpRfT2_Si2F, self).__init__(mesh, hybrid, 'inner', numbering_parameters, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_standard_inner_2_form')
        self._freeze_self_()




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_2/inner/main.py
    from __init__ import rfT2

    mesh = rfT2.rm(500, N_range=(3,3))

    f = mpRfT2_Si2F(mesh)


    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar
    import numpy as np
    def p(t, x, y): return np.sin(np.pi * x) * np.cos(np.pi * y) + t
    s = mpRfT2_Scalar(mesh, p)

    f.analytic_expression = s
    s.current_time = 0

    f.discretization()

    # R = f.reconstruct()

    f.visualization(show_mesh=True)
    # #
    print(f.error.L())
