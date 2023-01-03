# -*- coding: utf-8 -*-
"""
3 x 3 tensor field.


@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path:
    sys.path.append('./')

import numpy as np
from objects.CSCG._3d.fields.base import _3dCSCG_Continuous_FORM_BASE
from types import FunctionType, MethodType
from components.functions.timePlus3dSpace.constant import CFG
from functools import partial
from objects.CSCG._3d.fields.tensor.do.main import _3dCSCG_TensorField_DO
from objects.CSCG._3d.fields.tensor.numerical.main import _3dCSCG_TensorField_Numerical
from objects.CSCG._3d.fields.tensor.helpers.add import ___TENSOR_ADD_HELPER_1___
from objects.CSCG._3d.fields.tensor.helpers.sub import ___TENSOR_SUB_HELPER_1___
from objects.CSCG._3d.fields.tensor.helpers.neg import ___TENSOR_NEG_HELPER_1___

from objects.CSCG._3d.fields.tensor.visualize.main import _3dCSCG_TensorField_Visualize


class _3dCSCG_TensorField(_3dCSCG_Continuous_FORM_BASE, ndim=3):
    """The continuous vector field."""
    def __init__(self, mesh, func, ftype=None, valid_time=None, name='tensor-field'):
        if ftype is None:
            if isinstance(func, dict):
                ftype = 'boundary-wise'
            else:
                ftype = 'standard'
        else:
            pass
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_tensor_field')
        self.standard_properties.___PRIVATE_add_tag___('tensor_field')
        self.standard_properties.name = name

        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._DO_ = _3dCSCG_TensorField_DO(self)
        self._visualize_ = _3dCSCG_TensorField_Visualize(self)
        self._numerical_ = None
        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"3dCSCG_tensor_field=<{self.standard_properties.name}>@{id(self)}"

    def ___PRIVATE_set_func___(self, func, ftype='standard'):
        """
        Use this method to set up the function body and function type.

        Whenever define a new funcType, edit the currentFunc for the new type.
        """
        if ftype == 'standard':  # 3 by 3 functions in list or tuple.
            # standard func is function or method.
            assert np.shape(func) == (3, 3), f"Standard tensor only accepts list or tuple of shape (3,3)."

            _func_checked_ = [[None for _ in range(3)] for _ in range(3)]
            for i, fci_ in enumerate(func):
                for j, fij in enumerate(fci_):
                    if isinstance(fij, FunctionType):
                        # noinspection PyUnresolvedReferences
                        assert fij.__code__.co_argcount >= 4
                    elif isinstance(fij, MethodType):
                        # noinspection PyUnresolvedReferences
                        assert fij.__code__.co_argcount >= 5
                    elif isinstance(fij, (int, float)):
                        fij = CFG(fij)()

                    elif callable(fij):  # any other callable objects, we do not do check anymore.
                        pass

                    else:
                        raise Exception(f"standard tensor component [{i}][{j}] is a {fij.__class__}, wrong!")

                    # noinspection PyTypeChecker
                    _func_checked_[i][j] = fij

            self._func_ = _func_checked_

        else:
            raise Exception(f" <_3dCSCG_TensorField> does not accept funcType={ftype}")

        self._ftype_ = ftype

    def ___DO_evaluate_func_at_time___(self, time=None):
        """
        Evaluate the function at a particular time; reduce the number of variables from 4 to 3.

        :param float time: The time function is evaluated at.
        :return: A list of shape (3,) which can be sent to, for example, the instant function component of a form.
            They should be callable with ``(x,y,z)`` coordinates.
        :rtype: list
        """
        assert self.func is not None, f'Please first give me a func.'

        if time is None:
            time = self.current_time
        else:
            self.current_time = time

        if self._previous_func_id_time_[0:2] == (id(self.func), time):
            return self._previous_func_id_time_[2]
        else:
            if self.ftype == 'standard':
                RETURN_0 = partial(self.func[0][0], time), \
                    partial(self.func[0][1], time), \
                    partial(self.func[0][2], time)
                RETURN_1 = partial(self.func[1][0], time), \
                    partial(self.func[1][1], time), \
                    partial(self.func[1][2], time)
                RETURN_2 = partial(self.func[2][0], time), \
                    partial(self.func[2][1], time), \
                    partial(self.func[2][2], time)
                RETURN = [RETURN_0, RETURN_1, RETURN_2]

            else:
                raise Exception(f" do not understand funcType={self.ftype}")

            self._previous_func_id_time_ = (id(self.func), time, RETURN)

            return RETURN

    def reconstruct(self, *args, **kwargs):
        return self.do.reconstruct(*args, **kwargs)

    @property
    def shape(self):
        return 3, 3

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
            self._numerical_ = _3dCSCG_TensorField_Numerical(self)
        return self._numerical_

    def __neg__(self):
        """-self."""
        if self.ftype == 'standard':
            w0, w1, w2 = self.func
            w00, w01, w02 = w0
            w10, w11, w12 = w1
            w20, w21, w22 = w2

            x00 = ___TENSOR_NEG_HELPER_1___(w00)
            x01 = ___TENSOR_NEG_HELPER_1___(w01)
            x02 = ___TENSOR_NEG_HELPER_1___(w02)
            x10 = ___TENSOR_NEG_HELPER_1___(w10)
            x11 = ___TENSOR_NEG_HELPER_1___(w11)
            x12 = ___TENSOR_NEG_HELPER_1___(w12)
            x20 = ___TENSOR_NEG_HELPER_1___(w20)
            x21 = ___TENSOR_NEG_HELPER_1___(w21)
            x22 = ___TENSOR_NEG_HELPER_1___(w22)

            neg_tensor = _3dCSCG_TensorField(
                self.mesh,
                ([x00, x01, x02],
                 [x10, x11, x12],
                 [x20, x21, x22]),
                ftype='standard',
                valid_time=self.valid_time,
                name='-' + self.standard_properties.name
            )
            return neg_tensor

        else:
            raise Exception(f"cannot do neg for {self.ftype} _3dCSCG_TensorField.")

    def __sub__(self, other):
        """self - other"""
        if other.__class__.__name__ == '_3dCSCG_TensorField':

            if self.ftype == 'standard' and other.ftype == 'standard':

                w0, w1, w2 = self.func
                w00, w01, w02 = w0
                w10, w11, w12 = w1
                w20, w21, w22 = w2

                u0, u1, u2 = other.func
                u00, u01, u02 = u0
                u10, u11, u12 = u1
                u20, u21, u22 = u2

                x00 = ___TENSOR_SUB_HELPER_1___(w00, u00)
                x01 = ___TENSOR_SUB_HELPER_1___(w01, u01)
                x02 = ___TENSOR_SUB_HELPER_1___(w02, u02)
                x10 = ___TENSOR_SUB_HELPER_1___(w10, u10)
                x11 = ___TENSOR_SUB_HELPER_1___(w11, u11)
                x12 = ___TENSOR_SUB_HELPER_1___(w12, u12)
                x20 = ___TENSOR_SUB_HELPER_1___(w20, u20)
                x21 = ___TENSOR_SUB_HELPER_1___(w21, u21)
                x22 = ___TENSOR_SUB_HELPER_1___(w22, u22)

                sub_tensor = _3dCSCG_TensorField(
                    self.mesh,
                    ([x00, x01, x02],
                     [x10, x11, x12],
                     [x20, x21, x22]),
                    ftype='standard',
                    valid_time=self.valid_time,
                    name=self.standard_properties.name + '-' + other.standard_properties.name
                )
                return sub_tensor

            else:
                raise Exception(f"cannot do {self.ftype} _3dCSCG_TensorField - {other.ftype} _3dCSCG_TensorField")
        else:
            raise Exception(f"cannot do _3dCSCG_TensorField - {other.__class__}")

    def __add__(self, other):
        """self + other"""
        if other.__class__.__name__ == '_3dCSCG_TensorField':

            if self.ftype == 'standard' and other.ftype == 'standard':

                w0, w1, w2 = self.func
                w00, w01, w02 = w0
                w10, w11, w12 = w1
                w20, w21, w22 = w2

                u0, u1, u2 = other.func
                u00, u01, u02 = u0
                u10, u11, u12 = u1
                u20, u21, u22 = u2

                x00 = ___TENSOR_ADD_HELPER_1___(w00, u00)
                x01 = ___TENSOR_ADD_HELPER_1___(w01, u01)
                x02 = ___TENSOR_ADD_HELPER_1___(w02, u02)
                x10 = ___TENSOR_ADD_HELPER_1___(w10, u10)
                x11 = ___TENSOR_ADD_HELPER_1___(w11, u11)
                x12 = ___TENSOR_ADD_HELPER_1___(w12, u12)
                x20 = ___TENSOR_ADD_HELPER_1___(w20, u20)
                x21 = ___TENSOR_ADD_HELPER_1___(w21, u21)
                x22 = ___TENSOR_ADD_HELPER_1___(w22, u22)

                add_tensor = _3dCSCG_TensorField(
                    self.mesh,
                    ([x00, x01, x02],
                     [x10, x11, x12],
                     [x20, x21, x22]),
                    ftype='standard',
                    valid_time=self.valid_time,
                    name=self.standard_properties.name + '+' + other.standard_properties.name
                )
                return add_tensor

            else:
                raise Exception(f"cannot do {self.ftype} _3dCSCG_TensorField - {other.ftype} _3dCSCG_TensorField")
        else:
            raise Exception(f"cannot do _3dCSCG_TensorField + {other.__class__}")


if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\fields\tensor\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([1, 1, 2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto', 1), ('Lobatto', 1), ('Lobatto', 1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.cos(3*np.pi*z)

    TS = FC('tensor', ([p, 0, 0], [0, p, 0], [0, 0, p]))

    print(TS.shape)
