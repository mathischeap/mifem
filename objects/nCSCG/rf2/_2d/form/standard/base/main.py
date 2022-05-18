# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.form.base import _2nCSCG_RF2_FormBase
from objects.nCSCG.rf2._2d.form.standard.base.cochain import _2nCSCG_RF2_StandardFormCochain
from objects.nCSCG.rf2._2d.form.standard.base.IS import _2nCSCG_RF2_StandardFormIs


class _2nCSCG_RF2_StandardFormBase(_2nCSCG_RF2_FormBase):
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
        super(_2nCSCG_RF2_StandardFormBase, self).__init__(mesh, name)
        self.standard_properties.___PRIVATE_add_tag___('_2nCSCG_RF2_standard_form')
        assert orientation in ('inner', 'outer')
        self._orientation_ = orientation

        assert hybrid in (True,)
        self._hybrid_ = hybrid

        self._cochain_ = None
        self._IS_ = None

    #-------- must have methods ------------------------------------------------
    def ___Pr_check_func___(self, func):
        raise NotImplementedError()

    @property
    def orientation(self):
        return self._orientation_

    @property
    def cochain(self):
        if self._cochain_ is None:
            self._cochain_ = _2nCSCG_RF2_StandardFormCochain(self)
        return self._cochain_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _2nCSCG_RF2_StandardFormIs(self)
        return self._IS_







if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/form/standard/base/main.py
    pass
