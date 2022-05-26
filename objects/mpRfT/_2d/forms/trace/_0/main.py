# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:15 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.trace.base.main import mpRfT2_TraceFormBase
from objects.mpRfT._2d.forms.trace._0.numbering.main import mpRfT2_T0F_Numbering
from objects.mpRfT._2d.forms.trace._0.discretize.main import mpRfT2_T0F_Discretize
from objects.mpRfT._2d.forms.trace._0.reconstruct.main import mpRfT2_T0F_Reconstruct


class mpRfT2_T0F(mpRfT2_TraceFormBase):
    """"""
    def __init__(self, mesh, numbering_parameters='Naive', name='0-trace-form'):
        """

        Parameters
        ----------
        mesh
        name
        """
        super(mpRfT2_T0F, self).__init__(mesh, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_trace_0_form')

        self._numbering_ = mpRfT2_T0F_Numbering(self, numbering_parameters)
        self._discretize_ = mpRfT2_T0F_Discretize(self)
        self._reconstruct_ = mpRfT2_T0F_Reconstruct(self)

        self._freeze_self_()







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/trace/_0/main.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t0 = fc('0-t')

    print(t0.reconstruct)
