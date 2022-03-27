# -*- coding: utf-8 -*-
"""
Space related unittests.
"""
import sys
if './' not in sys.path: sys.path.append('/')
from root.config.main import *
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
import random

def test_Space_NO1_basis_functions_mapping_test():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print("sss {test_Space_NO1_basis_functions_mapping_test} ...... ", flush=True)

    if rAnk == mAster_rank:
        el1 = random.randint(1,3)
        el2 = random.randint(1,3)
        el3 = random.randint(1,3)
        c = random.uniform(0.0, 0.3)
        if c < 0.1:c = 0
    else:
        el1, el2, el3, c = [None for _ in range(4)]
    el1, el2, el3, c = cOmm.bcast([el1, el2, el3, c], root=mAster_rank)
    mesh = MeshGenerator('crazy', c=c)([el1, el2, el3])

    i, j, k = el3, el2, el1
    space = SpaceInvoker('polynomials')([('Lobatto', i), ('Lobatto', j), ('Lobatto', k)])
    FC = FormCaller(mesh, space)
    f2 = FC('2-f', is_hybrid=False)
    t2 = FC('2-t')

    ii, jj, kk = random.randint(4, 5), random.randint(4,6), random.randint(4,7)
    ratio = random.random()
    if ratio < 0.2:
        ratio = 0.2
    elif ratio > 0.9:
        ratio = 1
    else:
        pass
    xi, et, sg = np.linspace(-1,1,ii)*ratio, \
                 np.linspace(-1,1,jj)*ratio, \
                 np.linspace(-1,1,kk)*ratio

    co_N, bf2_N = f2.do.evaluate_basis_at_meshgrid([-1, ], et, sg, compute_xietasigma=True)
    bf2_N = bf2_N[0]
    bf2_N = bf2_N.reshape((i+1, j, k, jj*kk), order='F')
    bf2_N = bf2_N[0 ,...]
    bf2_N = bf2_N.reshape((j * k, jj * kk), order='F')

    co_S, bf2_S = f2.do.evaluate_basis_at_meshgrid([1, ], et, sg, compute_xietasigma=True)
    bf2_S = bf2_S[0]
    bf2_S = bf2_S.reshape((i+1, j, k, jj*kk), order='F')
    bf2_S = bf2_S[-1 ,...]
    bf2_S = bf2_S.reshape((j * k, jj * kk), order='F')

    co_W, bf2_W = f2.do.evaluate_basis_at_meshgrid(xi, [-1, ], sg, compute_xietasigma=True)
    bf2_W = bf2_W[1]
    bf2_W = bf2_W.reshape((i, j+1, k, ii*kk), order='F')
    bf2_W = bf2_W[:,0 ,...]
    bf2_W = bf2_W.reshape((i * k, ii * kk), order='F')

    co_E, bf2_E = f2.do.evaluate_basis_at_meshgrid(xi, [1, ], sg, compute_xietasigma=True)
    bf2_E = bf2_E[1]
    bf2_E = bf2_E.reshape((i, j+1, k, ii*kk), order='F')
    bf2_E = bf2_E[:,-1 ,...]
    bf2_E = bf2_E.reshape((i * k, ii * kk), order='F')

    co_B, bf2_B = f2.do.evaluate_basis_at_meshgrid(xi, et, [-1, ], compute_xietasigma=True)
    bf2_B = bf2_B[2]
    bf2_B = bf2_B.reshape((i, j, k+1, ii*jj), order='F')
    bf2_B = bf2_B[:,:,0 ,:]
    bf2_B = bf2_B.reshape((i * j, ii * jj), order='F')

    co_F, bf2_F = f2.do.evaluate_basis_at_meshgrid(xi, et, [1, ], compute_xietasigma=True)
    bf2_F = bf2_F[2]
    bf2_F = bf2_F.reshape((i, j, k+1, ii*jj), order='F')
    bf2_F = bf2_F[:,:,-1,:]
    bf2_F = bf2_F.reshape((i * j, ii * jj), order='F')

    co_T, bt2 = t2.do.evaluate_basis_at_meshgrid(xi, et, sg, compute_xietasigma=True)

    bt2_N = bt2['N'][0]
    bt2_S = bt2['S'][0]
    bt2_W = bt2['W'][0]
    bt2_E = bt2['E'][0]
    bt2_B = bt2['B'][0]
    bt2_F = bt2['F'][0]

    np.testing.assert_almost_equal(bf2_N, bt2_N)
    np.testing.assert_almost_equal(bf2_S, bt2_S)
    np.testing.assert_almost_equal(bf2_W, bt2_W)
    np.testing.assert_almost_equal(bf2_E, bt2_E)
    np.testing.assert_almost_equal(bf2_B, bt2_B)
    np.testing.assert_almost_equal(bf2_F, bt2_F)

    return 1

if __name__ == '__main__':
    # mpiexec -n 8 python _3dCSCG\TESTS\unittest_spaces.py
    test_Space_NO1_basis_functions_mapping_test()



