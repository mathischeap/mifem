# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 6:25 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from importlib import import_module
from objects.mpRfT._2d.forms.standard._1.base.numbering.do.main import mpRfT2_S1F_Numbering_DO
from objects.mpRfT._2d.forms.standard._1.base.numbering.local import mpRfT2_S1F_Numbering_Local


class mpRfT2_S1F_Numbering(FrozenOnly):
    """"""
    def __init__(self, f, numbering_parameters):
        """"""
        self._f_ = f
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
        self._numberer_ = getattr(import_module(path), scheme_name)(f)
        self._numbering_parameters_ = {'scheme_name': scheme_name}
        self._numbering_parameters_.update(parameters)
        self._gathering_ = None
        self._num_local_dofs_ = None
        self._do_ = mpRfT2_S1F_Numbering_DO(self)
        self._local_ = mpRfT2_S1F_Numbering_Local(self)
        self._freeze_self_()

    @property
    def local(self):
        return self._local_

    @property
    def gathering(self):
        if self._gathering_ is None:
            self._gathering_, self._num_local_dofs_ = self._numberer_(self._numbering_parameters_)
        return self._gathering_

    @property
    def GLOBAL_dofs(self):
        return self.gathering.GLOBAL_num_dofs

    @property
    def do(self):
        return self._do_



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
