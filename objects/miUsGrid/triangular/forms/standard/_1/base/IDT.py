# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/30 4:58 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

import numpy as np
from scipy.sparse import dia_matrix



class miUs_Triangular_S1F_InterfaceDofTopology(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._transitionMatrices_ = None
        self._transitionTypes_ = None
        self._freeze_self_()

    @property
    def EEs(self):
        raise NotImplementedError()

    @property
    def transition_matrices(self):
        """"""
        if self._transitionMatrices_ is None:

            EEs = self.EEs

            num_basis = self._sf_.num.basis

            tT = dict()
            tM = dict()

            key_pool = dict()

            for e in EEs:  # the elements that have transition matrices.

                key = ''.join([str(_) for _ in EEs[e]])

                if key in key_pool:
                    DIA = key_pool[key]

                else:

                    ones = np.ones(num_basis)
                    for edge_index in EEs[e]:
                        dofs = self._sf_.numbering.do.find.local_dofs_on_element_edge(edge_index)
                        ones[dofs] = -1

                    DIA = dia_matrix((ones, [0, ]), shape=(num_basis, num_basis))

                    key_pool[key] = DIA

                tM[e] = DIA
                tT[e] = key

            self._transitionMatrices_ = tM
            self._transitionTypes_ = tT

        return self._transitionMatrices_

    @property
    def transition_types(self):
        if self._transitionTypes_ is None:

            _ = self.transition_matrices

        return self._transitionTypes_

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
