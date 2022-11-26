# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 3:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.forms.standard._0.base.main import mpRfT2_S0F


class mpRfT2_Si0F(mpRfT2_S0F):
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
        super(mpRfT2_Si0F, self).__init__(mesh, hybrid, 'inner', numbering_parameters, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_standard_inner_0_form')
        self._freeze_self_()













if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_0/inner/main.py
    from __init__ import rfT2

    mesh = rfT2.mesh('crazy', c=0., bounds=([-2,2], [-2,2]))([15,15], 3)

    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar
    import numpy as np
    def p(t, x, y): return 1 / np.exp(np.abs(x**2 + y**2 - 2)) + t
    s = mpRfT2_Scalar(mesh, p)

    f = mpRfT2_Si0F(mesh)
    f.TW.func = s
    s.current_time = 0
    f.discretization()
    f.visualize()
    # mesh.visualize()

