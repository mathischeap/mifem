
import sys
if './' not in sys.path: sys.path.append('./')



from root.config.main import *
import random
from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_ColumnVector
from tools.linear_algebra.elementwise_cache.operators.concatenate.main import concatenate
from scipy import sparse as spspa


def test_LinearAlgebra_No0_EWC_ColumnVector():
    """"""
    if rAnk == mAster_rank:
        print("&&& [test_LinearAlgebra_No0_EWC_ColumnVector] ...... ", flush=True)

    if rAnk == mAster_rank:
        c = random.uniform(0, 0.15)
        i = random.randint(2, 3)
        j = random.randint(1, 3)
        k = random.randint(2, 4)
        l = random.randint(2, 3)
        m = random.randint(1, 3)
        n = random.randint(2, 4)
    else:
        c, i, j, k, l, m, n = None, None, None, None, None, None, None
    c, i, j, k, l, m, n = cOmm.bcast([c, i, j, k, l, m, n], root=mAster_rank)
    if c < 0.1: c = 0
    mesh = MeshGenerator('crazy', c=c)([i, j, k])
    space = SpaceInvoker('polynomials')([('Lobatto', l), ('Lobatto', m), ('Lobatto', n)])
    FC = FormCaller(mesh, space)

    def u(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def v(t, x, y, z): return np.sin(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2
    def w(t, x, y, z): return np.sin(np.pi*x) + np.cos(np.pi*y) * np.cos(np.pi*z-0.125)**2 + t
    def p(t, x, y, z): return np.cos(np.pi*x) + np.sin(np.pi*y) * np.sin(np.pi*z-0.125)**2 + t/2

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))
    f0 = FC('0-f', is_hybrid=False)
    f1 = FC('1-f', is_hybrid=False)
    f0.TW.func.body = scalar
    f0.TW.do.push_all_to_instant(0)
    f0.discretize()
    f1.TW.func.body = vector
    f1.TW.do.push_all_to_instant(0)
    f1.discretize()

    v0 = f0.cochain.EWC
    v1 = f1.cochain.EWC
    v2 = EWC_ColumnVector(mesh, 100)
    V = concatenate([v0, v1, v2])

    for i in V:
        v0_ = v0[i]
        v1_ = v1[i]
        v2_ = v2[i]
        V_ = V[i]
        vvv = spspa.vstack([v0_, v1_, v2_], format='csc')
        ZERO = V_ - vvv
        assert ZERO.nnz == 0

    return 1


if __name__ == '__main__':
    # mpiexec -n 10 python TESTS\unittest_linear_algebra.py

    test_LinearAlgebra_No0_EWC_ColumnVector()
