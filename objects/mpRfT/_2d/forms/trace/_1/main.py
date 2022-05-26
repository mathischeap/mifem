# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:15 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.trace.base.main import mpRfT2_TraceFormBase
from objects.mpRfT._2d.forms.trace._1.numbering.main import mpRfT2_T1F_Numbering
from objects.mpRfT._2d.forms.trace._1.discretize.main import mpRfT2_T1F_Discretize
from objects.mpRfT._2d.forms.trace._1.reconstruct.main import mpRfT2_T1F_Reconstruct


class mpRfT2_T1F(mpRfT2_TraceFormBase):
    """"""
    def __init__(self, mesh, numbering_parameters='Naive', name='1-trace-form'):
        """

        Parameters
        ----------
        mesh
        name
        """
        super(mpRfT2_T1F, self).__init__(mesh, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_trace_1_form')

        self._numbering_ = mpRfT2_T1F_Numbering(self, numbering_parameters)
        self._discretize_ = mpRfT2_T1F_Discretize(self)
        self._reconstruct_ = mpRfT2_T1F_Reconstruct(self)

        self._freeze_self_()





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/trace/_1/main.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t1 = fc('1-t')

    print(t1.cochain)
