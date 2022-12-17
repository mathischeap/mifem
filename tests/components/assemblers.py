# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/2/2022 4:31 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from root.config.main import RANK, MASTER_RANK
from components.assemblers import VectorAssembler
from components.distributors import VectorDistributor
import numpy as np


def test_vector_assembler():
    """"""
    if RANK == MASTER_RANK:
        print("ft2dw [test_vector_assembler] ...... ", flush=True)

    ADD = np.array([2., 3., 3., 3., 2., 1., 1., 1., 1., 1.])
    REP = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])

    G = {1: [0, 1, 2],
         2: [1, 2, 3],
         3: [4, 5, 6, 0, 3, 1, 2],
         '4': np.array([3, 4, 7, 8, 9])}
    D = {
        1: [1, 1, 1],
        2: [1, 1, 1],
        3: [1, 1, 1, 1, 1, 1, 1],
        '4': [1, 1, 1, 1, 1]
    }

    VA = VectorAssembler(G)
    output = VA(D, 'add')
    np.testing.assert_array_equal(output, ADD)
    output = VA(D, 'replace')
    np.testing.assert_array_equal(output, REP)

    VD = VectorDistributor(G)
    d = VD(output)
    for i in D:
        np.testing.assert_array_equal(d[i], D[i])

    G = [[0, 1, 2],
         [1, 2, 3],
         [4, 5, 6, 0, 3, 1, 2],
         np.array([3, 4, 7, 8, 9])]
    D = (
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1]
    )
    VA = VectorAssembler(G)
    output = VA(D, 'add')
    np.testing.assert_array_equal(output, ADD)
    output = VA(D, 'replace')
    np.testing.assert_array_equal(output, REP)

    VD = VectorDistributor(G)
    d = VD(output)
    for i, Di in enumerate(D):
        np.testing.assert_array_equal(d[i], Di)

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python tests/components/assemblers.py
    test_vector_assembler()
