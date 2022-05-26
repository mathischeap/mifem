# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/25/2022 9:18 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from importlib import import_module


class mpRfT2_T1F_Numbering(FrozenOnly):
    """"""
    def __init__(self, t, numbering_parameters):
        """"""
        self._t_ = t
        if isinstance(numbering_parameters, str):
            scheme_name = numbering_parameters
            parameters = dict()
        elif isinstance(numbering_parameters, dict): # do not use .pop() here
            scheme_name = numbering_parameters['scheme_name']
            parameters = dict()
            for key in numbering_parameters:
                if key != 'scheme_name':
                    parameters[key] = numbering_parameters[key]
        else:
            raise NotImplementedError()
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        path = base_path + scheme_name
        self._numberer_ = getattr(import_module(path), scheme_name)(t)
        self._numbering_parameters_ = {'scheme_name': scheme_name}
        self._numbering_parameters_.update(parameters)
        self.___Pr_reset_cache___()
        self._freeze_self_()

    def ___Pr_reset_cache___(self):
        """"""



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
