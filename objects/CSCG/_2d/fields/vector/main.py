# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

import numpy as np
from types import FunctionType, MethodType
# from SCREWS.frozen import FrozenOnly
# from BASE.elementwise_cache import EWC_ColumnVector
from objects.CSCG._2d.fields.base import _2dCSCG_Continuous_FORM_BASE
from functools import partial
from components.functions.timePlus2dSpace._0_ import _0t_
from components.functions.timePlus2dSpace.constant import CFGt as CFG_t_plus_2d
# from scipy import sparse as spspa
from objects.CSCG._2d.fields.vector.do.main import _2dCSCG_VectorField_DO
from objects.CSCG._2d.fields.vector.numerical.main import _2dCSCG_VectorField_Numerical
from objects.CSCG._2d.fields.vector.visualize.main import _2dCSCG_VectorField_Visualize

from objects.CSCG._2d.fields.scalar.helpers.mul import _2dCSCG_ScaMulHelper
from objects.CSCG._2d.fields.scalar.helpers.add import _2dCSCG_ScalarFieldAddHelper
from objects.CSCG._2d.fields.scalar.helpers.sub import _2dCSCG_ScalarFieldSubHelper
from objects.CSCG._2d.fields.scalar.helpers.neg import _2dCSCG_ScalarNeg


class _2dCSCG_VectorField(_2dCSCG_Continuous_FORM_BASE, ndim=2):
    """The continuous vector field."""
    def __init__(self, mesh, func, ftype=None, valid_time=None, name='vector-field'):
        """

        Parameters
        ----------
        mesh
        func
        ftype
        valid_time
        name
        """
        #--- when ftype is None, parse it from func type ----------------------------------------1
        if ftype is None:
            if isinstance(func, (list, tuple)):
                ftype = 'standard'
            elif isinstance(func, dict):
                ftype = 'boundary-wise'
            else:
                raise Exception()
        else:
            pass

        #----------------------------------------------------------------------------------------1
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_vector_field')
        self.standard_properties.___PRIVATE_add_tag___('vector_field')
        self.standard_properties.name = name
        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._do_ = _2dCSCG_VectorField_DO(self)
        self._numerical_ = _2dCSCG_VectorField_Numerical(self)
        self._visualize_ = _2dCSCG_VectorField_Visualize(self)
        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"2dCSCG_vector_field=<{self.standard_properties.name}>@{id(self)}"

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

        elif ftype == 'boundary-wise': # only valid (still as a vector) on mesh boundary (not domain boundary-wise)
            # no need to cover all mesh boundaries.
            assert isinstance(func, dict), f"when ftype == 'boundary-wise', " \
                                           f"we must put functions in a dict whose " \
                                           f"keys are boundary names and values are" \
                                           f"the functions."
            for bn in func:
                assert bn in self.mesh.boundaries.names, \
                    f"func key: [{bn}] is not a valid boundary name " \
                    f"({self.mesh.boundaries.names})"

                func_bn = func[bn]
                assert len(func_bn) == 2, \
                    f"2d vector should be of shape (2,), now it is {np.shape(func_bn)}."

                _func_bn_ck_ = list()
                for fci in func_bn:
                    # standard func is function or method.
                    if isinstance(fci, FunctionType):
                        assert fci.__code__.co_argcount >= 3
                    elif isinstance(fci, MethodType):
                        # noinspection PyUnresolvedReferences
                        assert fci.__code__.co_argcount >= 4
                    elif isinstance(fci, (int, float)):
                        fci = CFG_t_plus_2d(fci)()
                    else:
                        raise Exception()
                    _func_bn_ck_.append(fci)

                func[bn] = _func_bn_ck_

            self._func_ = func

        else:
            raise Exception(f" <_2dCSCG_VectorField> do not accept funcType={ftype}")
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

            elif self.ftype  == 'boundary-wise': #----------------------------------------------2
                RETURN = dict()
                for bn in self.func:
                    RETURN[bn] = [partial(self.func[bn][0], time),
                                  partial(self.func[bn][1], time)]

            else: #=============================================================================2
                raise Exception(f" do not understand funcType={self.ftype}")

            self._previous_func_id_time_ = (id(self.func), time, RETURN)
            return RETURN

    @property
    def shape(self):
        # noinspection PyRedundantParentheses
        return (2, )

    def reconstruct(self, *args, **kwargs):
        return self.do.reconstruct(*args, **kwargs)

    @property
    def do(self):
        return self._do_

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
                w0, w1 = self.func

                x0 = _2dCSCG_ScaMulHelper(w0, other)
                x1 = _2dCSCG_ScaMulHelper(w1, other)

                mul_vector = self.__class__(self.mesh,
                     [x0, x1],
                     ftype='standard',
                     valid_time=self.valid_time,
                     name=self.standard_properties.name + f'*{other}'
                )
                return mul_vector

            else:
                raise NotImplementedError()

        elif other.__class__.__name__ == '_2dCSCG_ScalarField':
            return other.__rmul__(self)

        else:
            raise NotImplementedError()

    def __rmul__(self, other):
        """other * self"""
        if other.__class__.__name__ in ('int', 'float', 'int64', 'int32'):
            if self.ftype == 'standard':
                w0, w1 = self.func

                x0 = _2dCSCG_ScaMulHelper(w0, other)
                x1 = _2dCSCG_ScaMulHelper(w1, other)

                mul_vector = self.__class__(self.mesh,
                     [x0, x1],
                     ftype='standard',
                     valid_time=self.valid_time,
                     name=f'{other}*' + self.standard_properties.name
                )
                return mul_vector

            else:
                raise NotImplementedError()

        elif other.__class__.__name__ == '_2dCSCG_ScalarField':
            return other.__mul__(self)

        else:
            raise NotImplementedError()



    def __add__(self, other):
        """self + other"""
        if other.__class__.__name__ == '_2dCSCG_VectorField':
            if self.ftype == 'standard' and other.ftype == 'standard':

                sf0, sf1 = self.func
                of0, of1 = other.func

                x0 = _2dCSCG_ScalarFieldAddHelper(sf0, of0)
                x1 = _2dCSCG_ScalarFieldAddHelper(sf1, of1)

                scalar = _2dCSCG_VectorField(self.mesh,
                     [x0, x1],
                     ftype='standard',
                     valid_time=self.valid_time,
                     name= self.standard_properties.name + '+' + other.standard_properties.name
                )
                return scalar

            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError(f"a 2d CSCG vector field can not + a {other.__class__.__name__}")

    def __sub__(self, other):
        """self - other"""
        if other.__class__.__name__ == '_2dCSCG_VectorField':
            if self.ftype == 'standard' and other.ftype == 'standard':

                sf0, sf1 = self.func
                of0, of1 = other.func

                x0 = _2dCSCG_ScalarFieldSubHelper(sf0, of0)
                x1 = _2dCSCG_ScalarFieldSubHelper(sf1, of1)

                scalar = _2dCSCG_VectorField(self.mesh,
                     [x0, x1],
                     ftype='standard',
                     valid_time=self.valid_time,
                     name= self.standard_properties.name + '-' + other.standard_properties.name
                )
                return scalar

            else:
                raise NotImplementedError()

        else:
            raise NotImplementedError(f"a 2d CSCG vector field can not - a {other.__class__.__name__}")


    def __neg__(self):
        """- self """
        if self.ftype == 'standard':

            x0 = _2dCSCG_ScalarNeg(self.func[0])
            x1 = _2dCSCG_ScalarNeg(self.func[1])

            scalar = _2dCSCG_VectorField(self.mesh,
                 [x0, x1],
                 ftype='standard',
                 valid_time=self.valid_time,
                 name= '-' + self.standard_properties.name
            )
            return scalar

        else:
            raise NotImplementedError()


if __name__ == '__main__':
    # mpiexec -n 6 python _2dCSCG\fields\vector\main.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * (x+1) * (y+1)
    def q(t, x, y): return t + np.cos(2*np.pi*x) * np.cos(np.pi*y) * (x+1) * (y+1)
    VV = FC('vector', [p,q])

    VV.visualize(time=1)
