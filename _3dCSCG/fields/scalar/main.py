# -*- coding: utf-8 -*-
"""
Continuous standard 3-form.

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

from types import FunctionType, MethodType
from _3dCSCG.fields.base import _3dCSCG_Continuous_FORM_BASE
from functools import partial
from root.config.main import *
from screws.functions.time_plus_3d_space.constant import CFG

from _3dCSCG.fields.scalar.do.main import _3dCSCG_ScalarField_DO
from _3dCSCG.fields.scalar.numerical import _3dCSCG_ScalarField_Numerical

from _3dCSCG.fields.scalar.helpers.neg import ___SCALAR_NEG_HELPER_1___
from _3dCSCG.fields.scalar.helpers.add import ___SCALAR_ADD_HELPER_1___
from _3dCSCG.fields.scalar.helpers.sub import ___SCALAR_SUB_HELPER_1___

from _3dCSCG.fields.scalar.visualize.main import _3dCSCG_ScalarField_Visualize


class _3dCSCG_ScalarField(_3dCSCG_Continuous_FORM_BASE, ndim=3):
    """The continuous scalar field."""
    def __init__(self, mesh, func, ftype=None, valid_time=None, name='scalar-field'):
        if ftype is None:
            if isinstance(func, dict):
                ftype= 'boundary-wise'
            else:
                ftype = 'standard'
        else:
            pass
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_scalar_field')
        self.standard_properties.name = name
        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._DO_ = _3dCSCG_ScalarField_DO(self)
        self._visualize_ = _3dCSCG_ScalarField_Visualize(self)
        self._numerical_ = None
        self._freeze_self_()

    def ___PRIVATE_set_func___(self, func, ftype='standard'):
        """
        Use this method to set up the function body and function type.

        Whenever define a new funcType, edit the currentFunc for the new type.
        """
        if ftype == 'standard':

            if isinstance(func, FunctionType):
                assert func.__code__.co_argcount >= 4
            elif isinstance(func, MethodType):
                # noinspection PyUnresolvedReferences
                assert func.__code__.co_argcount >= 5
            elif isinstance(func, (int, float)):
                func = CFG(func)()
            elif callable(func): # any other callable objects, we do not do check anymore.
                pass
            else:
                raise Exception()

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
                    assert func_bn.__code__.co_argcount >= 4
                elif isinstance(func_bn, MethodType):
                    # noinspection PyUnresolvedReferences
                    assert func_bn.__code__.co_argcount >= 5
                elif isinstance(func_bn, (int, float)):
                    func_bn = CFG(func_bn)()
                else:

                    raise Exception()

                self._func_[bn] = [func_bn,] # we always put a func representing a scalar in a list or tuple of shape (1,)
        elif ftype == 'trace-element-wise':
            # we have received a dict whose keys are local trace elements, values are callable that returns, xyz and a vector.
            assert isinstance(func, dict), f"func for trace-element-wise vector must a dict."
            for i in func: # valid local trace elements
                assert i in self.mesh.trace.elements, f"trace element #{i} is not in this core (#{rAnk})."
                # NOTE that we do not put the vector in a list or tuple, it should take (t, xi, eta, sigma) and then return xyz and the vector.
                assert callable(func[i]), f"func[{i}] is not callable."
            self._func_ = func
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

            elif self.ftype  == 'trace-element-wise':
                RETURN = dict()
                for i in self.func: # go through all valid trace elements
                    vi = self.func[i]
                    RETURN[i] = partial(vi, time) # We can see that for each trace-element, it is a single function

            else:
                raise Exception(f" Do not understand funcType={self.ftype}")

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
    def visualize(self):
        return self._visualize_
    @property
    def numerical(self):
        """The numerical property: A wrapper of all numerical methods, properties."""
        if self._numerical_ is None:
            self._numerical_ = _3dCSCG_ScalarField_Numerical(self)
        return self._numerical_

    def __neg__(self):
        """-self."""
        if self.ftype == 'standard':
            w0 = self.func[0]

            x0 = ___SCALAR_NEG_HELPER_1___(w0)

            neg_vector = _3dCSCG_ScalarField(self.mesh,
                                             x0,
                                             ftype='standard',
                                             valid_time=self.valid_time,
                                             name = '-' + self.standard_properties.name
                                            )
            return neg_vector

        else:
            raise Exception(f"cannot do neg for {self.ftype} _3dCSCG_ScalarField.")

    def __sub__(self, other):
        """self - other"""
        if other.__class__.__name__ == '_3dCSCG_ScalarField':

            if self.ftype == 'standard' and other.ftype == 'standard':

                w0 = self.func[0]
                u0 = other.func[0]

                x0 = ___SCALAR_SUB_HELPER_1___(w0, u0)

                sub_vector = _3dCSCG_ScalarField(self.mesh,
                                                 x0,
                                                 ftype='standard',
                                                 valid_time=self.valid_time,
                                                 name = self.standard_properties.name + '-' + other.standard_properties.name
                                                )
                return sub_vector

            else:
                raise Exception(f"cannot do {self.ftype} _3dCSCG_ScalarField - {other.ftype} _3dCSCG_ScalarField")
        else:
            raise Exception(f"cannot do _3dCSCG_ScalarField - {other.__class__}")

    def __add__(self, other):
        """self + other"""
        if other.__class__.__name__ == '_3dCSCG_ScalarField':

            if self.ftype == 'standard' and other.ftype == 'standard':

                w0 = self.func[0]
                u0 = other.func[0]

                x0 = ___SCALAR_ADD_HELPER_1___(w0, u0)

                add_vector = _3dCSCG_ScalarField(self.mesh,
                                                 x0,
                                                 ftype='standard',
                                                 valid_time=self.valid_time,
                                                 name = self.standard_properties.name + '+' + other.standard_properties.name
                                                )
                return add_vector

            else:
                raise Exception(f"cannot do {self.ftype} _3dCSCG_ScalarField - {other.ftype} _3dCSCG_ScalarField")
        else:
            raise Exception(f"cannot do _3dCSCG_ScalarField + {other.__class__}")











if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\fields\scalar\main.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([1,1,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.cos(3*np.pi*z)
    SS = FC('scalar', p)
    BS = FC('scalar', {'North': p, 'West':p})


    GV = SS.numerical.gradient

    BS.current_time=0
    BS.visualize()
