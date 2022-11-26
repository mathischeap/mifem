# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 12:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.forms.standard._1.base.main import mpRfT2_S1F
from objects.mpRfT._2d.forms.standard._1.inner.discretize.main import mpRfT2_Si1F_Discretize
from objects.mpRfT._2d.forms.standard._1.inner.migrate import mpRfT2_Si1F_Migrate





class mpRfT2_Si1F(mpRfT2_S1F):
    """"""

    def __init__(self, mesh, hybrid=True,
        numbering_parameters='Naive',  name='inner-oriented-1-form'):
        """

        Parameters
        ----------
        mesh
        hybrid :
            {True,}
        numbering_parameters
        name
        """
        super(mpRfT2_Si1F, self).__init__(mesh, hybrid, 'inner', numbering_parameters, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_standard_inner_1_form')

        self._discretize_ = mpRfT2_Si1F_Discretize(self)
        self._migrate_ = mpRfT2_Si1F_Migrate(self)
        self._freeze_self_()









if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/inner/main.py
    pass
