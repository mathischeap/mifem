# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 3:00 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from importlib import import_module


class mpRfT2_S1F_Coboundary(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    @property
    def incidence_matrix(self):
        return self._f_.matrices.incidence

    def __call__(self):
        """"""
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-4]) + '.'
        next_path = base_path + f'_2.{self._f_.orientation}.main'
        next_name = self._f_.__class__.__name__.replace('1', '2')
        next_class = getattr(import_module(next_path), next_name)
        next_form = next_class(self._f_.mesh, self._f_._hybrid_,
                               numbering_parameters=self._f_.numbering._numbering_parameters_,
                               name='d_'  + self._f_.standard_properties.name)
        LLC = self._f_.cochain.local
        E = self._f_.matrices.incidence
        dLLC = dict()
        for rp in LLC:
            dLLC[rp] = E[rp] @ LLC[rp]
        next_form.cochain.local = dLLC
        return next_form


if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/base/coboundary.py
    pass
