# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

import sys
if './' not in sys.path: sys.path.append('/')

from types import FunctionType, MethodType
from objects.CSCG._2d.fields.base import _2dCSCG_Continuous_FORM_BASE
from functools import partial
import numpy as np
from screws.functions.timePlus2dSpace._0_ import _0t_
from screws.functions.timePlus2dSpace.constant import CFGt as CFG_t_plus_2d
from objects.CSCG._2d.fields.scalar.do.main import _2dCSCG_ScalarField_DO
from objects.CSCG._2d.fields.scalar.numerical.main import _2dCSCG_ScalarField_Numerical
from objects.CSCG._2d.fields.scalar.visualize.main import _2dCSCG_ScalarField_Visualize


from objects.CSCG._2d.fields.scalar.helpers.mul import _2dCSCG_ScaMulHelper1, _2dCSCG_ScaMulHelper
from objects.CSCG._2d.fields.scalar.helpers.neg import _2dCSCG_ScalarNeg
from objects.CSCG._2d.fields.scalar.helpers.sub import _2dCSCG_ScalarFieldSubHelper
from objects.CSCG._2d.fields.scalar.helpers.add import _2dCSCG_ScalarFieldAddHelper


class _2dCSCG_ScalarField(_2dCSCG_Continuous_FORM_BASE, ndim=2):
    """The continuous scalar field."""
    def __init__(self, mesh, func, ftype=None, valid_time=None, name='scalar-field'):
        #--- when ftype is None, parse it from func type ----------------------------------------1
        if ftype is None:
            if callable(func) or isinstance(func, (int, float)):
                ftype = 'standard'
            elif isinstance(func, dict):
                ftype = 'boundary-wise'
            else:
                raise Exception()
        else:
            pass

        #----------------------------------------------------------------------------------------1
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_scalar_field')
        self.standard_properties.name = name
        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._DO_ = _2dCSCG_ScalarField_DO(self)
        self._numerical_ = _2dCSCG_ScalarField_Numerical(self)
        self._visualize_ = _2dCSCG_ScalarField_Visualize(self)
        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"2dCSCG_scalar_field=<{self.standard_properties.name}>@{id(self)}"

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

        elif ftype == 'boundary-wise': # mesh boundary wise (not domain boundary-wise)
            assert isinstance(func, dict), f" when ftype == 'boundary-wise', " \
                                           f"we must put functions in a dict whose " \
                                           f"keys are boundary names and values are" \
                                           f"the functions."

            self._func_ = dict()
            for bn in func:
                assert bn in self.mesh.boundaries.names, \
                    f"func key: [{bn}] is not a valid boundary name " \
                    f"({self.mesh.boundaries.names})"

                func_bn = func[bn]

                # standard func is function or method or int or float.
                if isinstance(func_bn, FunctionType):
                    assert func_bn.__code__.co_argcount >= 3
                elif isinstance(func_bn, MethodType):
                    # noinspection PyUnresolvedReferences
                    assert func_bn.__code__.co_argcount >= 4
                elif isinstance(func_bn, (int, float)):
                    func_bn = CFG_t_plus_2d(func_bn)()
                else:

                    raise Exception()

                self._func_[bn] = [func_bn,]
                # we always put a func representing a scalar in a list or tuple of shape (1,)

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

            elif self.ftype  == 'boundary-wise':

                RETURN = dict()
                for bn in self.func:
                    RETURN[bn] = [partial(self.func[bn][0], time),]

            else:
                raise Exception(f"Cannot do it for funcType={self.ftype}")

            self._previous_func_id_time_ = (id(self.func), time, RETURN)
            return RETURN

    @property
    def shape(self):
        return 1,

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


    def __mul__(self, other):
        """self * other"""
        if other.__class__.__name__ in ('int', 'float', 'int64', 'int32'):
            if self.ftype == 'standard':
                w0 = self.func[0]

                x0 = _2dCSCG_ScaMulHelper(w0, other)

                mul_vector = _2dCSCG_ScalarField(self.mesh,
                         x0,
                         ftype='standard',
                         valid_time=self.valid_time,
                         name=self.standard_properties.name + f'*{other}'
                         )
                return mul_vector
            else:
                raise NotImplementedError()

        elif other.__class__.__name__ == '_2dCSCG_VectorField':
            if self.ftype == 'standard' and other.ftype == 'standard':

                sfunc = self.func[0]
                vf0, vf1 = other.func

                x0 = _2dCSCG_ScaMulHelper1(sfunc, vf0)
                x1 = _2dCSCG_ScaMulHelper1(sfunc, vf1)

                mul_vector = other.__class__(self.mesh,
                     [x0, x1],
                     ftype='standard',
                     valid_time=self.valid_time,
                     name=self.standard_properties.name + '*' + other.standard_properties.name
                     )
                return mul_vector

            else:
                raise NotImplementedError()

        else:
            raise NotImplementedError()

    def __rmul__(self, other):
        """other * self"""
        if other.__class__.__name__ in ('int', 'float', 'int64', 'int32'):
            if self.ftype == 'standard':
                w0 = self.func[0]

                x0 = _2dCSCG_ScaMulHelper(w0, other)

                mul_vector = _2dCSCG_ScalarField(self.mesh,
                                 x0,
                                 ftype='standard',
                                 valid_time=self.valid_time,
                                 name= f'{other}*' + self.standard_properties.name
                                 )
                return mul_vector
            else:
                raise NotImplementedError()

        elif other.__class__.__name__ == '_2dCSCG_VectorField':
            if self.ftype == 'standard' and other.ftype == 'standard':

                sfunc = self.func[0]
                vf0, vf1 = other.func

                x0 = _2dCSCG_ScaMulHelper1(sfunc, vf0)
                x1 = _2dCSCG_ScaMulHelper1(sfunc, vf1)

                mul_vector = other.__class__(self.mesh,
                     [x0, x1],
                     ftype='standard',
                     valid_time=self.valid_time,
                     name=other.standard_properties.name + '*' + self.standard_properties.name
                     )
                return mul_vector

            else:
                raise NotImplementedError()

        else:
            raise NotImplementedError()


    def __add__(self, other):
        """self + other"""
        if other.__class__.__name__ == '_2dCSCG_ScalarField':
            if self.ftype == 'standard' and other.ftype == 'standard':

                sf = self.func[0]
                of = other.func[0]

                x0 = _2dCSCG_ScalarFieldAddHelper(sf, of)

                scalar = _2dCSCG_ScalarField(self.mesh,
                                 x0,
                                 ftype='standard',
                                 valid_time=self.valid_time,
                                 name= self.standard_properties.name + '+' + other.standard_properties.name
                                 )
                return scalar

            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError(f"a 2d CSCG scalar field can not + a {other.__class__.__name__}")

    def __sub__(self, other):
        """self - other"""
        if other.__class__.__name__ == '_2dCSCG_ScalarField':
            if self.ftype == 'standard' and other.ftype == 'standard':

                sf = self.func[0]
                of = other.func[0]

                x0 = _2dCSCG_ScalarFieldSubHelper(sf, of)

                scalar = _2dCSCG_ScalarField(self.mesh,
                                 x0,
                                 ftype='standard',
                                 valid_time=self.valid_time,
                                 name= self.standard_properties.name + '-' + other.standard_properties.name
                                 )
                return scalar

            else:
                raise NotImplementedError()

        else:
            raise NotImplementedError(f"a 2d CSCG scalar field can not - a {other.__class__.__name__}")


    def __neg__(self):
        """- self """
        if self.ftype == 'standard':

            x0 = _2dCSCG_ScalarNeg(self.func[0])

            scalar = _2dCSCG_ScalarField(self.mesh,
                             x0,
                             ftype='standard',
                             valid_time=self.valid_time,
                             name= '-' + self.standard_properties.name
                             )
            return scalar

        else:
            raise NotImplementedError()





if __name__ == '__main__':
    # mpiexec -n 6 python _2dCSCG\field\scalar.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y)
    SS = FC('scalar', p)
    SS.current_time = 100

    SS.visualize()

