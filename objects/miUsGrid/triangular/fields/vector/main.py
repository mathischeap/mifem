# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 1:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.fields.base import miUsGrid_TriangularFieldBase

from functools import partial
from types import FunctionType, MethodType
from screws.functions.time_plus_2d_space._0_ import _0t_

from objects.miUsGrid.triangular.fields.vector.reconstruct.main import miUsGrid_Triangular_Vector_Reconstruct
from objects.miUsGrid.triangular.fields.vector.do.main import miUsGrid_Triangular_Vector_Do
from objects.miUsGrid.triangular.fields.vector.numerical.main import miUsGrid_Triangular_Vector_Numerical


class miUsGrid_Triangular_Vector(miUsGrid_TriangularFieldBase):
    """"""

    def __init__(self, mesh, func, ftype=None, valid_time=None, name='Tri-Vector'):
        """"""
        super(miUsGrid_Triangular_Vector, self).__init__(mesh, valid_time, name)

        #--- when ftype is None, parse it from func type ----------------------------------------1
        if ftype is None:
            if isinstance(func, (list, tuple)):
                ftype = 'standard'
            else:
                raise Exception()
        else:
            pass

        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)

        self._reconstruct_ = miUsGrid_Triangular_Vector_Reconstruct(self)
        self._do_ = miUsGrid_Triangular_Vector_Do(self)
        self._numerical_ = miUsGrid_Triangular_Vector_Numerical(self)

        self._freeze_self_()


    def __repr__(self):
        """"""
        return f"miUsGrid_Triangular_Vector_Field=<{self.standard_properties.name}>@{id(self)}"


    def ___PRIVATE_set_func___(self, func, ftype='standard'):
        """Use this method to set up the function body and function type.

        Whenever define a new funcType, edit the currentFunc for the new type.
        """
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
            raise Exception(f" <miUsGrid_Triangular_Vector> do not accept funcType={ftype}")
        self._ftype_ = ftype

    def ___DO_evaluate_func_at_time___(self, time=None):
        """Evaluate the function at a particular time; reduce the number of variables from 4 to 3.

        :param float time: The time function is evaluated at.
        :return: A list of shape (3,) which can be sent to, for example, the instant function
            component of a form. They should be callable with `(x,y,z)` coordinates.
        :rtype: list

        """
        if time is None:
            time = self.current_time
        else:
            self.current_time = time

        self.current_time = time
        assert self.func is not None, 'Please first set func.'
        if self._previous_func_id_time_[0:2] == (id(self.func), time):
            return self._previous_func_id_time_[2]
        else:
            if self.ftype == 'standard': #------------------------------------------------------2
                RETURN = partial(self.func[0], time), partial(self.func[1], time)

            else: #=============================================================================2
                raise Exception(f" do not understand funcType={self.ftype}")

            self._previous_func_id_time_ = (id(self.func), time, RETURN)
            return RETURN

    @property
    def shape(self):
        return 2,

    @property
    def numerical(self):
        return self._numerical_



if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/fields/vector/main.py
    def fx(t, x, y): return t + x - y
    def fy(t, x, y): return t - x + y
    from objects.miUsGrid.triangular.__test__.Random.test_mesh import mesh

    v = miUsGrid_Triangular_Vector(mesh, (fx, fy))
    v.current_time = 0

    import numpy as np
    xi= et = np.linspace(-1,1,10)
    R = v.reconstruct(xi, et)
