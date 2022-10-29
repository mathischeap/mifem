# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/9 14:28
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


class nLS_Customize(FrozenOnly):
    def __init__(self, nLS):
        """"""
        self._nLS_ = nLS
        self.___customizations___ = list()
        self._freeze_self_()

    @property
    def customizations(self):
        """The current existing customizations."""
        return self.___customizations___

    def clear(self, i=None):
        """Clear the ith customization. If `i` is None, clear all.

        Parameters
        ----------
        i

        Returns
        -------

        """
        if i is None:
            self.___customizations___ = dict()
        else:
            raise NotImplementedError(f"clear i={i} customization is not implemented.")

    def set_no_evaluation(self, i):
        """Let the nonlinear system do not affect the value of #r dof.

        So dof_i will be equal to dof_0 (the initial value (or initial guess)).
        """
        self.___customizations___.append(('set_no_evaluation', i))

    def set_fixed_solution(self, i, value):
        """ So we make the initial value of dof #`i` to be `value`, plus we do `self.set_no_evaluation(i)`.

        Parameters
        ----------
        i :
            The dof #i.
        value :
            The value of dof #'i' will be this one.

        Returns
        -------

        """
        self.___customizations___.append(('set_ith_value_of_initial_guess', [i, value]))
        self.___customizations___.append(('set_no_evaluation', i))

    def apply_strong_BC(self, i, j, dof_itp, cochain_itp=None, AS='local'):
        """We will only need partial dofs, we do not need to set particular values other than 0
        for the right hand vector.

        Parameters
        ----------
        i
        j
        dof_itp
        cochain_itp
        AS

        Returns
        -------

        """
        if hasattr(dof_itp, 'standard_properties') and 'form' in dof_itp.standard_properties.tags:
            dof_itp = dof_itp.BC.interpret
        else:
            pass


        # check 1 _______________________________________________________________2
        if i == j:
            assert cochain_itp is None, f"when i == j is None, we must have pc is None."

        # check 2 _______________________________________________________________2
        if cochain_itp is None:
            assert i == j, \
                f"when do not provide pc, we must set diagonal block, " \
                f"so i == j, now i={i}, j={j}."
        elif hasattr(cochain_itp, 'standard_properties') and 'form' in cochain_itp.standard_properties.tags:
            cochain_itp = cochain_itp.BC.interpret
        else:
            pass

        #======== customize =============================================================
        I, J = self._nLS_.shape
        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        assert j % 1 == 0, f"j={j}({j.__class__.__name__}) cannot be an index!"
        assert 0 <= i < I and 0 <= j < J, f"(i,j)= ({i},{j}) is out of index range!"

        if i == j:
            self.___customizations___.append((
                'set_no_evaluations',       # key
                [i, dof_itp, AS]            # para
            ))
        else:
            raise NotImplementedError(dof_itp, cochain_itp)








if __name__ == '__main__':
    # mpiexec -n 4 python tools/linear_algebra/nonlinear_system/customize.py
    pass
