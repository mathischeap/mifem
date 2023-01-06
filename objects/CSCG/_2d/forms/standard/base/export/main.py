# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/15 3:46 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')

from components.freeze.base import FrozenOnly


class _2dCSCG_Standard_Form_Export(FrozenOnly):
    """Export the standard form to a file."""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        raise NotImplementedError()


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/base/export/main.py
    from __init__ import cscg2

    mesh = cscg2.mesh('crazy')([2, 2, 2])
