# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/12 4:20 PM
"""
import sys

import numpy as np

if './' not in sys.path: sys.path.append('./')
from screws.miscellaneous.miprint import miprint
import screws.miscellaneous.mirand as random
from tools.linear_algebra.gathering.irregular.ir_chain_matrix.main import iR_Chain_Gathering_Matrix
from __init__ import rfT2
from tools.linear_algebra.elementwise_cache.operators.bmat.main import bmat

def test_mpRfT2_ir_numbering():
    """"""
    miprint(f"-sr- [test_mpRfT2_ir_numbering] ... ", flush=True)
    fc = rfT2.rf(100)

    f2 = fc('2-f-o')
    f1 = fc('1-f-o')
    nt = fc('nst')

    gm2 = f2.numbering.gathering
    gm1 = f1.numbering.gathering
    gmt = nt.numbering.gathering

    GM = iR_Chain_Gathering_Matrix([gm2, gmt, gm1])
    num_dofs_f2 = f2.num.GLOBAL_dofs
    num_dofs_t = nt.num.GLOBAL_dofs
    assert not GM.___Pr_IS_regular___

    for rp in GM:
        GV = GM[rp]
        gv2 = gm2[rp]
        gvt = gmt[rp]
        gv1 = gm1[rp]
        gv = np.concatenate([gv2, gvt+num_dofs_f2, gv1+num_dofs_f2+num_dofs_t])
        np.testing.assert_array_equal(GV, gv)

    NDfs = GM.GLOBAL_num_dofs
    dofs = random.sample(range(NDfs), int(NDfs/100))
    for dof in dofs:
        EI = GM.do.find.elements_and_local_indices_of_dof(dof)
        if EI is None:
            pass
        else:
            for rp, i in zip(*EI):
                assert GM[rp][i] == dof, f"must be the case!"

    GM = iR_Chain_Gathering_Matrix(gm1)
    for rp in GM:
        GV = GM[rp]
        np.testing.assert_array_equal(GV, gm1[rp])

    M = f1.matrices.mass
    E = f1.matrices.incidence

    A = bmat([[M, E.T ],
              [E, None]])

    A.gathering_matrices = ([f1, f2], [f1, f2])

    A = A.assembled

    return 1


if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/__tests__/unittests/numbering.py
    test_mpRfT2_ir_numbering()
