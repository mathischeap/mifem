# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:21 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.form.standard.base.main import _2nCSCG_RF2_StandardFormBase
from objects.nCSCG.rf2._2d.form.standard._0.base.numbering.main import _2nCSCG_RF2_Standard0Form_Numbering
from objects.nCSCG.rf2._2d.form.standard._0.base.discretize.main import _2nCSCG_RF2_S0F_Discretize
from objects.nCSCG.rf2._2d.form.standard._0.base.reconstruct.main import _2nCSCG_RF2_S0F_Reconstruct
from objects.nCSCG.rf2._2d.form.standard._0.base.visualize import _2nCSCG_RF2_S0F_Visualize
from objects.nCSCG.rf2._2d.form.standard._0.base.error import _2nCSCG_RF2_S0F_Error
from objects.nCSCG.rf2._2d.form.standard._0.base.do import _2nCSCG_RF2_S0F_Do


class _2nCSCG_RF2_Standard0FormBase(_2nCSCG_RF2_StandardFormBase):
    """"""

    def __init__(self, mesh, hybrid, orientation, numbering_parameters, name):
        """"""
        super(_2nCSCG_RF2_Standard0FormBase, self).__init__(mesh, hybrid, orientation, name)
        self.standard_properties.___PRIVATE_add_tag___('_2nCSCG_RF2_standard_0_form')
        self._k_ = 0

        self._numbering_ = _2nCSCG_RF2_Standard0Form_Numbering(self, numbering_parameters)
        self._discretize_ = _2nCSCG_RF2_S0F_Discretize(self)
        self._reconstruct_ = _2nCSCG_RF2_S0F_Reconstruct(self)
        self._visualize_ = None
        self._error_ = None
        self._do_ = None

    #-------- must have methods ------------------------------------------------
    def ___Pr_check_func___(self, func):
        assert func.mesh is self.mesh

        if func.__class__.__name__ == '_2nCSCG_RF2_ScalarField':
            assert func.ftype in ('standard',), \
                f"2nCSCG RF2 0form FUNC do not accept func _2nCSCG_RF2_ScalarField of ftype {func.ftype}."
        else:
            raise Exception(f"2nCSCG RF2 0form FUNC do not accept func {func.__class__}")

    @property
    def numbering(self):
        return self._numbering_

    @property
    def discretize(self):
        return self._discretize_

    @property
    def reconstruct(self):
        return self._reconstruct_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2nCSCG_RF2_S0F_Visualize(self)
        return self._visualize_

    @property
    def error(self):
        if self._error_ is None:
            self._error_ = _2nCSCG_RF2_S0F_Error(self)
        return self._error_

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _2nCSCG_RF2_S0F_Do(self)
        return self._do_


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
