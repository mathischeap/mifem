# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/10 2:03 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from components.miscellaneous.miprint import miprint
from tools.linearAlgebra.gathering.chain import GatheringMatrixChaining
from __init__ import miTri, tools
import numpy as np
from root.config.main import RANK, MASTER_RANK

def _func(t, x, y): return np.sin(2*np.pi*x) * np.sin(2*np.pi*y) + t
def _fx(t, x, y): return np.sin(2*np.pi*x) * np.cos(np.pi*y) + t
def _fy(t, x, y): return np.cos(2*np.pi*x) * np.sin(np.pi*y) + t


class CustomizeSequent(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("--- [GM-Customize-Sequent-Test] ...... ", flush=True)
        self._freeze_self_()


    def __call__(self):
        """"""
        fc = miTri.call('rand0', 2)
        f0 = fc('0-f-i')
        f1o = fc('1-f-o')
        f1i = fc('1-f-i')
        f2 = fc('2-f-i')

        scalar = fc('scalar', _func)
        vector = fc('vector', (_fx, _fy))

        scalar.current_time = 0
        vector.current_time = 0

        f0.CF = scalar
        f1i.CF = vector
        f1o.CF = vector
        f2.CF = scalar

        f0.discretize()
        f1i.discretize()
        f1o.discretize()
        f2.discretize()
        M0 = f0.matrices.mass
        M1i = f1i.matrices.mass
        M1o = f1o.matrices.mass
        M2 = f2.matrices.mass

        CGM = GatheringMatrixChaining(f0, f1o, f1i, f2)(chain_method='sequent')

        A = ([M0  , None, None, None],
             [None, M1i , None, None],
             [None, None, M1o , None],
             [None, None, None, M2  ])

        A = tools.linalg.bmat(A)
        A.assembler.chain_method = 'sequent'
        A.gathering_matrices = (CGM, CGM)

        f0.BC.CF = scalar
        f0.BC.boundaries = ['Upper','Right']
        f1i.BC.CF = vector
        f1i.BC.boundaries = ['Down','Left']
        f1o.BC.CF = vector
        f1o.BC.boundaries = ['Right','Upper','Down']

        A.customize.identify_global_rows_according_to(0, f0.BC.interpret)
        A.customize.identify_global_rows_according_to(1, f1i.BC.interpret)
        A.customize.identify_global_rows_according_to(2, f1o.BC.interpret)

        A = A.assembled
        A = A.do.gather_M_to_core()

        if RANK == MASTER_RANK:
            A = A.tocsr()
            for i in range(A.shape[0]):
                Ai = A[i,:]
                nnz = Ai.nnz

                if nnz == 1:
                    index = Ai.indices[0]
                    assert i == index

        return 1



if __name__ == "__main__":
    # mpiexec -n 4 python tests/tools/gatheringMatrix/customize_sequent.py
    CustomizeSequent()()
