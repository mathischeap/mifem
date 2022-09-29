# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/27/2022 12:35 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from screws.miscellaneous.miprint import miprint
from screws.miscellaneous.mirand import randint

from __init__ import miTri


class miUsGrid_Triangle_StandardFormNumbering(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("TriSFN [miUsGrid_Triangle_StandardFormNumbering] ...... ", flush=True)
        p = randint(2, 8)
        self.fc = miTri.form('test0', p)
        self._freeze_self_()


    def __call__(self):
        return 1

if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/__test__/unittests/standard_forms/numbering.py
    miUsGrid_Triangle_StandardFormNumbering()()
