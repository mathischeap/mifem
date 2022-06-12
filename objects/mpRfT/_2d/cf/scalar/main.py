# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 3:51 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from functools import partial
from types import FunctionType, MethodType
from screws.functions.time_plus_2d_space._0_ import _0t_

from objects.mpRfT._2d.cf.base import mpRfT2_ContinuousField
from objects.mpRfT._2d.cf.scalar.reconstruct import mpRfT2_ScalarReconstruct
from objects.mpRfT._2d.cf.scalar.visualize import mpRfT2_ScalarVisualize





class mpRfT2_Scalar(mpRfT2_ContinuousField):
    """"""
    def __init__(self, mesh, func, ftype='standard', valid_time=None, name='scalar-field'):
        super().__init__(mesh, ftype, valid_time, name)
        self.___Pr_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._reconstruct_ = None
        self._visualize_ = None
        self._freeze_self_()

    #---------- built-in ---------------------------------------------------------------
    def ___Pr_set_func___(self, func, ftype='standard'):
        """"""
        if ftype == 'standard':
            if isinstance(func, FunctionType):
                # noinspection PyUnresolvedReferences
                assert func.__code__.co_argcount >= 3
            elif isinstance(func, MethodType):
                # noinspection PyUnresolvedReferences
                assert func.__code__.co_argcount >= 4
            elif isinstance(func, (int, float)) and func == 0:
                func = _0t_
            elif callable(func): # any other callable objects, we do not do check anymore.
                pass
            else:
                raise Exception(f"func={func} is a {func.__class__}, cannot be understood.")
            self._func_ = [func,]

        else:
            raise Exception(f" <ScalarField> do not accept funcType={ftype}")

        self._ftype_ = ftype

    def ___Pr_evaluate_func___(self):
        """
        Evaluate the function at a particular time; reduce the number of variables from 4 to 3.

        :return: A list of shape (1,) which can be sent to the instant function component of a form.
        :rtype: list
        """
        time = self.current_time
        assert self.func is not None, 'Please first set func.'

        if self._previous_func_id_time_[0:2] == (id(self.func), time):
            RETURN = self._previous_func_id_time_[2]

        else:
            if self.ftype == 'standard':
                RETURN = [partial(self.func[0], time),]
            else:
                raise Exception(f"Cannot do it for funcType={self.ftype}")

            self._previous_func_id_time_ = (id(self.func), time, RETURN)

        return RETURN

    @property
    def shape(self):
        return 1,

    @property
    def ___parameters___(self):
        return None

    #------------------------------- personal -------------------------------------------
    @property
    def reconstruction(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = mpRfT2_ScalarReconstruct(self)
        return self._reconstruct_

    @property
    def visualization(self):
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_ScalarVisualize(self)
        return self._visualize_







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/cf/scalar/main.py
    from __init__ import rfT2
    import numpy as np

    mesh = rfT2.rm(50)

    def p(t, x, y): return np.sin(2*np.pi*x) * np.cos(2 * np.pi*y)  + t

    s = mpRfT2_Scalar(mesh, p)
    s.current_time = 0

    s.visualization()