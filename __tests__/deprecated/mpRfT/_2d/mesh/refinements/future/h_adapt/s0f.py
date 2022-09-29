# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/22 1:28 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import cOmm, MPI, np
from screws.collections.single_list import single_list



class hAdapt_S0F(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, approach=1, **kwargs):
        """"""
        if approach == 1:
            self.___A1___(**kwargs)
        else:
            raise NotImplementedError(f"a_adapt_to s0f approach = {approach} is not implemented.")

    def ___A1___(self, levels):
        """"""
        #------- check levels -------------------------------------------------------------
        if isinstance(levels, (int, float)):
            levels = np.array([levels,])
        else:
            levels = list(levels)
            levels.sort()
            levels = np.array(levels)
        assert np.ndim(levels) == 1 and np.min(levels) > 0 and np.max(levels) < 1 \
               and np.all(np.diff(levels)) > 0, \
            f"levels={levels} wrong, must be 1d, and in (0,1), and increasing."
        #-------------------------------------------------------------------------------------------
        f = self._f_
        mesh = f.mesh
        num_levels = len(levels)
        coo = mesh.coo_map.uniform(2**num_levels * 5, ndim=1)
        V = f.reconstruction(coo, value_only=True)
        V = V.bcW # basic-cell-wise data struction.
        MIN = list()
        MAX = list()
        for i in V:
            Vi = np.abs(V[i])
            MIN.append(np.min(Vi))
            MAX.append(np.max(Vi))
        MIN = np.min(MIN)
        MAX = np.max(MAX)
        MIN = cOmm.allreduce(MIN, op=MPI.MIN)
        MAX = cOmm.allreduce(MAX, op=MPI.MAX)
        delta= MAX - MIN
        levels = [delta * _ for _ in levels]
        rdf = self.___Pr_hAdapt_from_bcW___(V, levels)
        mesh.refinements.future._rfd_ = rdf

    def ___Pr_hAdapt_from_bcW___(self, bcW, levels):
        """"""
        RDF = list()
        for i in bcW:

            V =bcW[i]
            ik = str(i) + '-'

            rdf = self.___Pr_divide_cell___(V, levels, ik)

            if isinstance(rdf, str):
                pass
            else:
                rdf = single_list(rdf)
                RDF.extend(rdf)

        RDF = dict.fromkeys(RDF, self._f_.mesh.dN)

        return RDF

    def ___Pr_divide_cell___(self, V, levels, ik):
        """"""
        if len(levels) == 0 or np.max(np.abs(V)) < levels[0]:
            return ik
        else:
            Vs = self.___Pr_divide_V_into_4_parts___(V)
            lv0_hR = any([np.average(_) > levels[0] for _ in Vs])
            if not lv0_hR:
                return ik
            else:
                rdf = [ik + str(_) for _ in range(4)]
                for j, vs in enumerate(Vs):
                    output = self.___Pr_divide_cell___(vs, levels[1:], rdf[j])
                    rdf[j] = output
                return rdf

    @staticmethod
    def ___Pr_divide_V_into_4_parts___(V):
        """"""
        S = V.shape[0]
        P = int(S / 2)
        return V[:P, :P], V[P:, :P], V[:P, P:], V[P:, P:]







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/refinements/future/h_adapt/s0f.py
    from __init__ import rfT2
    from objects.mpRfT._2d.forms.standard._0.inner.main import mpRfT2_Si0F
    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar

    mesh = rfT2.mesh('crazy', c=0.0, bounds=([-2, 2], [-2, 2]))([15, 15], 3)



    def p(t, x, y): return 1 / np.exp(np.abs(x**2 + y**2 - 2))**2 + 0 * t
    s = mpRfT2_Scalar(mesh, p)
    f = mpRfT2_Si0F(mesh)
    f.TW.func = s
    s.current_time = 0
    f.discretization()
    # f.visualization(show_mesh=True)
    mesh.refinements.future.h_adapt_to(f, levels=(0.5, 0.75))
    mesh = mesh.do.evolve()
    mesh.visualize()