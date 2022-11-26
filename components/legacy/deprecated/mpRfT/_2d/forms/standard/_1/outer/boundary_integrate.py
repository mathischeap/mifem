# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/01 5:47 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
import numpy as np
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_SparseMatrix
from scipy.sparse import csr_matrix


class mpRfT2_So1F_BI(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, t, **kwargs):
        """"""
        if t.__class__.__name__ == 'mpRfT2_NSgF':
            return self.___Pr_with_mpRfT2_NSgF___(t, **kwargs)
        else:
            raise NotImplementedError(f"Boundary integral for mpRfT2_So1F with "
                                      f"{t.__class__.__name__} is not implemented.")

    def ___Pr_with_mpRfT2_NSgF___(self, t, qdp=2):
        """

        Parameters
        ----------
        t
        qdp : int
            Quadrature degree plus.

        Returns
        -------

        """
        f = self._f_
        mesh = f.mesh

        fBasis = dict()
        tBasis = dict()
        TDict = dict()

        Qnw = mesh.space.Gauss

        BI = dict()

        for crp in mesh.rcfc:
            cell = mesh[crp]
            fN = f.N[crp]

            Mrc = list()


            if fN in TDict:
                TU, TD, TL, TR = TDict[fN]
            else:
                T = cell.space.trace_matrix._2dCSCG_1Trace_Outer[1]
                TU, TD, TL, TR = T['U'].T, T['D'].T, T['L'].T, T['R'].T
                TDict[fN] = [TU, TD, TL, TR]

            Frame = cell.frame
            for edge in Frame:
                segments = Frame[edge]

                for seg in segments:
                    tN = t.N[seg]
                    o, d = seg.origin_and_delta
                    if isinstance(o, tuple):
                        if seg.direction == 'UD':
                            o = o[0]
                        else:
                            o = o[1]
                    else:
                        pass

                    n, w = Qnw(tN + qdp)

                    kf = (tN, fN, o, d)
                    if kf in fBasis:
                        fb = fBasis[kf]
                    else:
                        n_sg = o + (n + 1) * d / 2
                        fb = mesh.space[fN].basises[0].edge_basis(x=n_sg)
                        fBasis[kf] = fb

                    if tN in tBasis:
                        tb = tBasis[tN]
                    else:
                        tb = mesh.space[tN].basises[0].node_basis(x=n)
                        tBasis[tN] = tb

                    M_seg = np.einsum('ik, jk, k -> ij', fb, tb, w, optimize='greedy')

                    if edge == 'U':
                        T = TU
                    elif edge == 'D':
                        T = TD
                    elif edge == 'L':
                        T = TL
                    else:
                        T = TR

                    M_seg = T @ M_seg
                    Mrc.append(M_seg)

            BI[crp] = csr_matrix(np.hstack(Mrc))

        return EWC_SparseMatrix(self._f_.mesh.rcfc,
                                BI,
                                cache_key_generator='no_cache'
                                )







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/outer/boundary_integrate.py
    from __init__ import rfT2

    fc = rfT2.rf(100, N_range=(3,3))

    f = fc('1-f-o')
    t = fc('nst')

    # import numpy as np
    # def p(t, x, y): return np.sin(np.pi*x) * np.cos(np.pi*y) + t
    # def q(t, x, y): return np.cos(np.pi*x) * np.sin(np.pi*y) + t
    #
    # def h(t, x, y): return np.sin(np.pi*x) * np.sin(np.pi*y) + t
    #
    # s = fc('scalar', h)
    # v = fc('vector', (p, q))

    # f.TW.func = v
    # v.current_time = 0
    # f.discretize()

    # t.TW.func = s
    # s.current_time = 0
    # t.discretize()
    # f.mesh.visualize()
    B = f.boundary_integrate(t)

    M = f.matrices.mass

    MB = M @ B

    print(MB)