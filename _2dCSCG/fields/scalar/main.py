# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

import sys
if './' not in sys.path: sys.path.append('./')

from types import FunctionType, MethodType
from _2dCSCG.fields.base import _2dCSCG_Continuous_FORM_BASE
from functools import partial
import numpy as np
from screws.functions.time_plus_2d_space._0_ import _0t_
from _2dCSCG.fields.scalar.do.main import _2dCSCG_ScalarField_DO
from _2dCSCG.fields.scalar.numerical.main import _2dCSCG_ScalarField_Numerical
from _2dCSCG.fields.scalar.visualize.main import _2dCSCG_ScalarField_Visualize






class _2dCSCG_ScalarField(_2dCSCG_Continuous_FORM_BASE, ndim=2):
    """The continuous scalar field."""
    def __init__(self, mesh, func, ftype='standard', valid_time=None, name='scalar-field'):
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_scalar_field')
        self.standard_properties.name = name
        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._DO_ = _2dCSCG_ScalarField_DO(self)
        self._numerical_ = _2dCSCG_ScalarField_Numerical(self)
        self._visualize_ = _2dCSCG_ScalarField_Visualize(self)
        self._freeze_self_()


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
            raise Exception(f" <ScalarField> do not accept funcType={ftype}")
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
        return (1,)

    def reconstruct(self, *args, **kwargs):
        return self.do.reconstruct(*args, **kwargs)

    @property
    def do(self):
        return self._DO_

    @property
    def numerical(self):
        return self._numerical_

    @property
    def visualize(self):
        return self._visualize_




if __name__ == '__main__':
    # mpiexec -n 6 python _2dCSCG\field\scalar.py
    from _2dCSCG.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y)
    SS = FC('scalar', p)
    SS.current_time = 100

    SS.visualize()

