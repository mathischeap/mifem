# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 7:43 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')


from functools import partial
from types import FunctionType, MethodType
from components.functions.timePlus2dSpace._0_ import _0t_

from objects.mpRfT._2d.cf.base import mpRfT2_ContinuousField
from objects.mpRfT._2d.cf.vector.reconstruct import mpRfT2_VectorReconstruct
from objects.mpRfT._2d.cf.vector.visualize import mpRfT2_VectorVisualize





class mpRfT2_Vector(mpRfT2_ContinuousField):
    """"""

    def __init__(self, mesh, func, ftype='standard', valid_time=None, name='vector-field'):
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
            # standard func is function or method.
            assert len(func) == 2
            func = list(func)
            for i, fci in enumerate(func):
                if isinstance(fci, FunctionType):
                    # noinspection PyUnresolvedReferences
                    assert fci.__code__.co_argcount >= 3
                elif isinstance(fci, MethodType):
                    # noinspection PyUnresolvedReferences
                    assert fci.__code__.co_argcount >= 4
                elif isinstance(fci, (int, float)) and fci == 0:
                    func[i] = _0t_
                elif callable(fci): # any other callable objects, we do not do check anymore.
                    pass
                else:
                    raise Exception()
            self._func_ = func
        else:
            raise Exception(f" <_2dCSCG_VectorField> do not accept funcType={ftype}")
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
            return self._previous_func_id_time_[2]
        else:
            if self.ftype == 'standard':
                RETURN = partial(self.func[0], time), partial(self.func[1], time)
            else:
                raise Exception(f" do not understand funcType={self.ftype}")
            self._previous_func_id_time_ = (id(self.func), time, RETURN)
            return RETURN

    @property
    def shape(self):
        return 2,

    @property
    def ___parameters___(self):
        return None

    #-------------- personal -----------------------------------------------------------
    @property
    def reconstruction(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = mpRfT2_VectorReconstruct(self)
        return self._reconstruct_

    @property
    def visualization(self):
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_VectorVisualize(self)
        return self._visualize_






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/cf/vector/main.py
    from __init__ import rfT2

    mesh = rfT2.rm(100)

    import numpy as np
    def p(t, x, y): return np.sin(np.pi*x) * np.sin(np.pi*y) + t
    def q(t, x, y): return np.cos(np.pi*x) * np.cos(np.pi*y) + t

    v = mpRfT2_Vector(mesh, (p,q))
    v.current_time = 0

    v.visualization()
