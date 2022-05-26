# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:12 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.base import mpRfT2_FormBase

from objects.mpRfT._2d.forms.trace.base.cochain import mpRfT2_TF_Cochain
from objects.mpRfT._2d.forms.trace.base.IS import mpRfT2_TF_IS


class mpRfT2_TraceFormBase(mpRfT2_FormBase):
    """"""

    def __init__(self, mesh, name):
        """

        Parameters
        ----------
        mesh
        name
        """
        super(mpRfT2_TraceFormBase, self).__init__(mesh, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_trace_form')

        self._IS_ = mpRfT2_TF_IS(self)
        self._cochain_ = mpRfT2_TF_Cochain(self)
        self._numbering_ = None
        self._num_ = None

        self._discretize_ = None
        self._reconstruct_ = None

        self._migrate_ = None


    #-------- must have methods ------------------------------------------------
    def ___Pr_check_func___(self, func):
        """"""
        assert func.mesh is self.mesh

        if func.__class__.__name__ == 'mpRfT2_Scalar':
            assert func.ftype in ('standard',), \
                f"mpRfT2-trace-form FUNC do not accept func mpRfT2_Scalar of ftype {func.ftype}."
        else:
            raise Exception(f"mpRfT2-trace-form FUNC do not accept func {func.__class__}")


    @property
    def cochain(self):
        return self._cochain_

    @property
    def IS(self):
        return self._IS_

    @property
    def numbering(self):
        return self._numbering_

    @property
    def num(self):
        return self._num_

    @property
    def discretize(self):
        return self._discretize_

    @property
    def reconstruct(self):
        return self._reconstruct_

    @property
    def migrate(self):
        return self._migrate_









if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
