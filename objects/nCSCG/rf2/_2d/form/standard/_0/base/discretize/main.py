# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.form.standard._0.base.discretize.scalar.standard import \
    _2nCSCG_RF2_S0F_Discretize_Standard_Scalar






class _2nCSCG_RF2_S0F_Discretize(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self.___Pr_dis_standard_scalar___ = _2nCSCG_RF2_S0F_Discretize_Standard_Scalar(f)
        self._freeze_self_()


    def __call__(self, update_cochain=True, target='func'):
        """

        Parameters
        ----------
        update_cochain : bool
            If we send the output to `cochain.local`.
        target : str

        Returns
        -------

        """

        if target == 'func':

            if self._f_.TW.func.__class__.__name__ == '_2nCSCG_RF2_ScalarField':

                if self._f_.TW.func.ftype == 'standard':
                    LCC =  self.___Pr_dis_standard_scalar___(target)
                else:
                    raise NotImplementedError(f"2nCSCG 0-form cannot (target func) "
                                              f"discretize _2nCSCG_RF2_ScalarField of ftype={self._f_.TW.func.ftype}")

            else:
                raise NotImplementedError(f'2nCSCG 0-form can not (target func) '
                                          f'discretize {self._f_.TW.func.__class__}.')

        else:
            raise NotImplementedError(f"2dCSCG outer 0-form cannot discretize "
                                      f"while targeting at {target}.")

        if update_cochain: self._f_.cochain.local = LCC

        return LCC




if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/form/standard/_0/base/discretize/main.py
    pass
