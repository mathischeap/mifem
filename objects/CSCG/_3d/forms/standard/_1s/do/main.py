# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/05 9:55 PM
"""
from objects.CSCG._3d.forms.standard.base.do import _3dCSCG_Standard_Form_DO


class _3dCSCG_S1F_Do(_3dCSCG_Standard_Form_DO):
    """"""

    def __init__(self, sf):
        """"""
        super(_3dCSCG_S1F_Do, self).__init__(sf)
        self._freeze_self_()
