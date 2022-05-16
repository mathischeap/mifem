# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 7:43 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')


from functools import partial
from types import FunctionType, MethodType
from screws.functions.time_plus_2d_space._0_ import _0t_

from objects.nCSCG.rf2._2d.fields.base import _2nCSCG_FieldBase
from objects.nCSCG.rf2._2d.fields.vector.reconstruct import _2nCSCG_RF2_VectorReconstruct
from objects.nCSCG.rf2._2d.fields.vector.visualize import _2nCSCG_RF2_VectorVisualize


class _2nCSCG_RF2_VectorField(_2nCSCG_FieldBase):
    """"""

    def __init__(self, mesh, func, ftype='standard', valid_time=None, name='vector-field'):
        super().__init__(mesh, ftype, valid_time, name)
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_RF2_Vector')
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

    #-------------- personal --------------------------------------------------------------------
    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = _2nCSCG_RF2_VectorReconstruct(self)
        return self._reconstruct_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2nCSCG_RF2_VectorVisualize(self)
        return self._visualize_



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
