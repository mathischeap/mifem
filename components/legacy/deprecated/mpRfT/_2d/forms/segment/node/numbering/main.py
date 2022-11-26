# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/25/2022 9:18 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

import numpy as np
from components.freeze.base import FrozenOnly
from importlib import import_module
from tools.linearAlgebra.gathering.vector import Gathering_Vector
from tools.linearAlgebra.gathering.irregular.ir_matrix.main import iR_Gathering_Matrix
from objects.mpRfT._2d.forms.segment.node.numbering.do.main import mpRfT2_NSgF_Numbering_DO
from objects.mpRfT._2d.forms.segment.node.numbering.local import mpRfT2_NSgF_Numbering_Local


class mpRfT2_NSgF_Numbering(FrozenOnly):
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
        self._numbering_parameters_ = {'scheme_name': scheme_name,}
        self._numbering_parameters_.update(parameters)
        self._gathering_ = None
        self._sgW_GM_ = None
        self._num_local_dofs_ = None
        self._do_ = mpRfT2_NSgF_Numbering_DO(self)
        self._local_ = mpRfT2_NSgF_Numbering_Local(self)
        self._freeze_self_()

    @property
    def local(self):
        return self._local_

    @property
    def gathering(self):
        if self._gathering_ is None:
            sgWGM = self.sgW_gathering
            mesh = self._t_.mesh

            GVs = dict()
            for rp in mesh.rcfc:
                cell = mesh[rp]
                Frame = cell.frame

                vector = list()

                for edge in Frame:
                    segments = Frame[edge]
                    for seg in segments:
                        srp = seg.__repr__()
                        vector.append(sgWGM[srp].full_vector)

                vector = np.concatenate(vector)
                GVs[rp] = Gathering_Vector(rp, vector)

            self._gathering_ = iR_Gathering_Matrix(GVs, mesh_type='mpRfT2')

        return self._gathering_

    @property
    def sgW_gathering(self):
        if self._sgW_GM_ is None:
            self._sgW_GM_, self._num_local_dofs_ = \
                self._numberer_(self._numbering_parameters_)
        return self._sgW_GM_

    @property
    def GLOBAL_dofs(self):
        return self.gathering.GLOBAL_num_dofs

    @property
    def do(self):
        return self._do_




if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/node/numbering/main.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t = fc('nst')

    print(t.num.local_dofs)
