# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/25/2022 9:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.gathering.vector import Gathering_Vector
from tools.linearAlgebra.gathering.irregular.ir_matrix.main import iR_Gathering_Matrix
from root.config.main import SIZE, RANK, COMM



class Naive(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._freeze_self_()


    def __call__(self, parameters):
        """"""
        scheme_name = parameters['scheme_name']
        assert scheme_name == 'Naive', f"parameters wrong."

        #--- no other parameters, we use the default numbering scheme ----------
        if len(parameters) == 1:
            rcWGM = self.___Pr_no_para_routine___()
        else:
            raise NotImplementedError(f"Scheme {parameters} not implemented.")

        return rcWGM

    def ___Pr_no_para_routine___(self):
        """A stupid routine."""
        mesh = self._t_.mesh
        NUMb = self._t_.num.basis

        if RANK == 0:
            numbered = dict()
            START = 0
            for seg in mesh.segments:
                # go through all local segments in RANK #0
                num_basis = NUMb[seg]
                numbered[seg.__repr__()] = range(START, START + num_basis)
                START += num_basis

        else:
            numbered, START = COMM.recv(source=RANK - 1, tag=RANK)
            for seg in mesh.segments:
                rp = seg.__repr__()
                if rp in numbered:
                    pass
                else:
                    num_basis = NUMb[seg]
                    numbered[rp] = range(START, START + num_basis)
                    START += num_basis

        if RANK != SIZE - 1:
            COMM.send([numbered, START], dest=RANK + 1, tag=RANK + 1)

        GVs = dict()
        NUM_local_dofs = 0
        for seg in mesh.segments:
            rp = seg.__repr__()
            GVs[rp] = Gathering_Vector(rp, numbered[rp])
            NUM_local_dofs += NUMb[seg]

        # noinspection PyUnboundLocalVariable
        rcWGM = iR_Gathering_Matrix(GVs, mesh_type='mpRfT2')
        # noinspection PyUnboundLocalVariable
        return rcWGM, NUM_local_dofs







if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/node/numbering/Naive.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t = fc('nst')

    print(t.numbering.sgW_gathering)
