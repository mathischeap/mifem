# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 12:36 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.standard.base.main import mpRfT2_StandardFormBase

from objects.mpRfT._2d.forms.standard._2.base.numbering.main import mpRfT2_S2F_Numbering
from objects.mpRfT._2d.forms.standard._2.base.num import mpRfT2_S2F_Num
from objects.mpRfT._2d.forms.standard._2.base.error import mpRfT2_S2F_Error


from objects.mpRfT._2d.forms.standard._2.base.discretize.main import mpRfT2_S2F_Discretize
from objects.mpRfT._2d.forms.standard._2.base.reconstruct import mpRfT2_S2F_Reconstruct
from objects.mpRfT._2d.forms.standard._2.base.migrate import mpRfT2_S2F_Migrate


from objects.mpRfT._2d.forms.standard._2.base.visaulize import mpRfT2_S2F_Visualize
from objects.mpRfT._2d.forms.standard._2.base.matrices.main import mpRfT2_S2F_Matrices


class mpRfT2_S2F(mpRfT2_StandardFormBase):
    """"""

    def __init__(self, mesh, hybrid, orientation, numbering_parameters, name):
        """"""
        super(mpRfT2_S2F, self).__init__(mesh, hybrid, orientation, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_standard_2_form')
        self._k_ = 2

        self._numbering_ = mpRfT2_S2F_Numbering(self, numbering_parameters)
        self._num_ = mpRfT2_S2F_Num(self)
        self._error_ = mpRfT2_S2F_Error(self)

        self._discretize_ = mpRfT2_S2F_Discretize(self)
        self._reconstruct_ = mpRfT2_S2F_Reconstruct(self)
        self._migrate_ = mpRfT2_S2F_Migrate(self)

        self._visualize_ = mpRfT2_S2F_Visualize(self)
        self._matrices_ = mpRfT2_S2F_Matrices(self)


    def ___Pr_check_func___(self, func):
        """"""
        assert func.mesh is self.mesh

        if func.__class__.__name__ == 'mpRfT2_Scalar':
            assert func.ftype in ('standard',), \
                f"mpRfT2_S2F FUNC do not accept func mpRfT2_Scalar of ftype {func.ftype}."
        else:
            raise Exception(f"mpRfT2_S2F FUNC do not accept func {func.__class__}")


    @property
    def visualize(self):
        return self._visualize_

    @property
    def matrices(self):
        return self._matrices_






if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
