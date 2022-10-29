# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 6:27 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.gathering.vector import Gathering_Vector
from tools.linearAlgebra.gathering.irregular.ir_matrix.main import iR_Gathering_Matrix
from root.config.main import SIZE, RANK, COMM, np


class Naive(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()


    def __call__(self, parameters):
        """"""
        scheme_name = parameters['scheme_name']
        assert scheme_name == 'Naive', f"parameters wrong"

        #--- no other parameters, we use the default numbering scheme ----------
        if len(parameters) == 1:
            rcWGM = self.___Pr_no_para_routine___()
        else:
            raise NotImplementedError(f"Scheme {parameters} not implemented.")

        return rcWGM

    def ___Pr_no_para_routine___(self):
        """"""
        assert self._f_.IS.hybrid, f"this routine only works for hybrid mpRfT2 standard-1-form."
        mesh = self._f_.mesh
        for rank in range(SIZE):
            if rank == RANK:
                if RANK == 0:
                    current = 0
                else:
                    current = COMM.recv(source=RANK - 1, tag=RANK)

                MY_START = current

                for rp in mesh.rcfc:
                    num_basis = self._f_.num.basis[rp]
                    current += num_basis

                if RANK < SIZE-1:
                    COMM.send(current, dest=RANK + 1, tag=RANK + 1)

                GVs = dict()
                NUM_local_dofs = 0
                for rp in mesh.rcfc:
                    num_basis = self._f_.num.basis[rp]
                    GVs[rp] = Gathering_Vector(rp, np.arange(MY_START, MY_START + num_basis))
                    MY_START += num_basis
                    NUM_local_dofs += num_basis

        # noinspection PyUnboundLocalVariable
        rcWGM = iR_Gathering_Matrix(GVs, mesh_type='mpRfT2')
        # noinspection PyUnboundLocalVariable
        return rcWGM, NUM_local_dofs

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
