# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/7 12:19
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class nLS_DO(FrozenOnly):
    def __init__(self, nLS):
        """"""
        self._nLS_ = nLS
        self._freeze_self_()

    def evaluate_f(self, unknown_variables, neg=False):
        """

        Parameters
        ----------
        unknown_variables : list
            Usually be the shadows of all unknown_variables.

        neg : bool


        Returns
        -------

        """
        S0, S1 = self._nLS_.shape
        f = list()

        #------ linear terms contribution ---------------------------------------------
        for i in range(S0):
            fi = list()
            for j in range(S1):
                Aij = self._nLS_.A[i][j]
                if Aij == 0 or Aij is None:
                    pass
                else:
                    # noinspection PyUnresolvedReferences
                    fi.append(Aij @ unknown_variables[j])
            f.append(fi)

        #-------- nonlinear terms contribution --------------------------------------------
        f_nt = self._nLS_.nonlinear_terms.___PRIVATE_evaluate___(unknown_variables)

        #---------- right hand side vector contribution -----------------------------------
        v = self._nLS_.b

        #============ combine them below ==================================================
        for i in range(S0):
            f[i].extend(f_nt[i])
            if v[i] is not None:
                f[i].append(-v[i])

            LENfi = len(f[i])

            if LENfi == 0:
                raise Exception(f"Cannot be. A trivial check.")
            elif LENfi == 1:
                # noinspection PyUnresolvedReferences
                f[i] = f[i][0]
            elif LENfi == 2:
                # noinspection PyUnresolvedReferences
                f[i] = f[i][0] + f[i][1]
            else:
                fi0 = f[i][0]
                f[i] = fi0.___PRIVATE_sum___(f[i][1:])

        if neg:
            for i in range(S0):
                f[i] = - f[i]
        else:
            pass

        return f



if __name__ == '__main__':
    # mpiexec -n 4 python tools/linear_algebra/nonlinear_system/do.py
    pass
