# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 3:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.forms.base import mpRfT2_FormBase
from objects.mpRfT._2d.forms.standard.base.IS import mpRfT2_SF_IS
from objects.mpRfT._2d.forms.standard.base.cochain import mpRfT2_SF_Cochain
from objects.mpRfT._2d.forms.standard.base.N import mpRfT2_SF_N


class mpRfT2_StandardFormBase(mpRfT2_FormBase):
    """"""
    def __init__(self, mesh, hybrid, orientation, name):
        """

        Parameters
        ----------
        mesh
        hybrid :
            {True,}
            - True: Fully hybrid, each sub-cell is separated.

        """
        super(mpRfT2_StandardFormBase, self).__init__(mesh, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_standard_form')
        assert orientation in ('inner', 'outer')
        self._orientation_ = orientation

        assert hybrid in (True,)
        self._hybrid_ = hybrid

        self._IS_ = mpRfT2_SF_IS(self)
        self._cochain_ = mpRfT2_SF_Cochain(self)
        self._N_ = mpRfT2_SF_N(self)

        self._numbering_ = None
        self._num_ = None
        self._error_ = None
        self._discretize_ = None
        self._reconstruct_ = None
        self._migrate_ = None
        self._visualize_ = None



    def ___Pr_check_analytic_expression___(self, func):
        raise NotImplementedError()

    def ___Pr_check_BC_analytic_expression___(self, ae):
        raise NotImplementedError()

    @property
    def orientation(self):
        return self._orientation_

    @property
    def IS(self):
        return self._IS_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def N(self):
        return self._N_



    @property
    def numbering(self):
        """"""
        return self._numbering_

    @property
    def num(self):
        """"""
        return self._num_

    @property
    def error(self):
        return self._error_

    @property
    def discretization(self):
        return self._discretize_

    @property
    def reconstruction(self):
        return self._reconstruct_

    @property
    def migration(self):
        return self._migrate_

    @property
    def visualization(self):
        return self._visualize_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
