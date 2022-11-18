# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 1:56 PM
"""
import sys

import numpy as np

if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.fields.base import miUsGrid_TriangularFieldBase

from functools import partial
from types import FunctionType, MethodType
from screws.functions.timePlus2dSpace._0_ import _0t_

from objects.miUsGrid.triangular.fields.scalar.reconstruct.main import miUsGrid_Triangular_Scalar_Reconstruct
from objects.miUsGrid.triangular.fields.scalar.do.main import miUsGrid_Triangular_Scalar_Do
from objects.miUsGrid.triangular.fields.scalar.numerical.main import miUsGrid_Triangular_Scalar_Numerical
from objects.miUsGrid.triangular.fields.scalar.export.main import miUsGrid_Triangular_Scalar_Export


class miUsGrid_Triangular_Scalar(miUsGrid_TriangularFieldBase):
    """"""

    def __init__(self, mesh, func, ftype=None, valid_time=None, name='Tri-Scalar'):
        """"""
        super(miUsGrid_Triangular_Scalar, self).__init__(mesh, valid_time, name)

        #--- when ftype is None, parse it from func type ----------------------------------------1
        if ftype is None:
            if callable(func) or isinstance(func, (int, float)):
                ftype = 'standard'
            else:
                raise Exception()
        else:
            pass
        #----------------------------------------------------------------------------------------1
        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)

        self._reconstruct_ = miUsGrid_Triangular_Scalar_Reconstruct(self)
        self._do_ = miUsGrid_Triangular_Scalar_Do(self)
        self._numerical_ = miUsGrid_Triangular_Scalar_Numerical(self)
        self._export_ = None

        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"miUsGrid_Triangular_Scalar_Field=<{self.standard_properties.name}>@{id(self)}"


    def ___PRIVATE_set_func___(self, func, ftype='standard'):
        """
        Use this method to set up the function body and function type.

        Whenever define a new funcType, edit the currentFunc for the new type.
        """
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
            raise Exception(f" <miUsGrid_Triangular_Scalar> do not accept funcType={ftype}")

        self._ftype_ = ftype

    def ___DO_evaluate_func_at_time___(self, time=None):
        """
        Evaluate the function at a particular time; reduce the number of variables from 4 to 3.

        :param float time: The time function is evaluated at.
        :return: A list of shape (1,) which can be sent to the instant function component of a form.
        :rtype: list
        """
        if time is None:
            time = self.current_time
        else:
            self.current_time = time

        assert self.func is not None, 'Please first set func.'
        if self._previous_func_id_time_[0:2] == (id(self.func), time):
            return self._previous_func_id_time_[2]

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
    def numerical(self):
        return self._numerical_

    @property
    def export(self):
        if self._export_ is None:
            self._export_ = miUsGrid_Triangular_Scalar_Export(self)
        return self._export_


if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/fields/scalar/main.py
    def func(t, x, y): return t + x + y
    from tests.objects.miUsGrid.triangular.randObj.rand_mesh import mesh

    s = miUsGrid_Triangular_Scalar(mesh, func)
    s.current_time = 0

    xi= et = np.linspace(-1,1,10)
    R = s.reconstruct(xi, et)