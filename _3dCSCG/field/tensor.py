# -*- coding: utf-8 -*-
"""
3 x 3 tensor field.


@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

import numpy as np
from SCREWS.frozen import FrozenOnly
from _3dCSCG.field.main import _3dCSCG_Continuous_FORM_BASE
from types import FunctionType, MethodType
from SCREWS.functions._4d import CFG
from functools import partial

from importlib import import_module
from SCREWS.numerical._4d import NumericalPartialDerivative_txyz_Functions


class _3dCSCG_TensorField(_3dCSCG_Continuous_FORM_BASE, ndim=3):
    """The continuous vector field."""
    def __init__(self, mesh, func, ftype=None, valid_time=None, name='tensor-field'):
        if ftype is None:
            if isinstance(func, dict):
                ftype= 'boundary-wise'
            else:
                ftype = 'standard'
        else:
            pass
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_tensor_field')
        self.standard_properties.name = name

        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._DO_ = _3dCSCG_TensorField_DO(self)
        self._numerical_ = None
        self._freeze_self_()


    def ___PRIVATE_set_func___(self, func, ftype='standard'):
        """
        Use this method to set up the function body and function type.

        Whenever define a new funcType, edit the currentFunc for the new type.
        """
        if ftype == 'standard': # 3 by 3 functions in list or tuple.
            # standard func is function or method.
            assert np.shape(func) == (3,3), f"Standard tensor only accepts list or tuple of shape (3,3)."
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
        if time is None:
            time = self.current_time
        else:
            self.current_time = time

        assert self.func is not None, 'Please first set func.'

        if self._previous_func_id_time_[0:2] == (id(self.func), time):
            return self._previous_func_id_time_[2]
        else:
            if self.ftype == 'standard':
                RETURN_0 = partial(self.func[0][0], time), partial(self.func[0][1], time), partial(self.func[0][2], time)
                RETURN_1 = partial(self.func[1][0], time), partial(self.func[1][1], time), partial(self.func[1][2], time)
                RETURN_2 = partial(self.func[2][0], time), partial(self.func[2][1], time), partial(self.func[2][2], time)
                RETURN = [RETURN_0, RETURN_1, RETURN_2]

            else:
                raise Exception(f" do not understand funcType={self.ftype}")

            self._previous_func_id_time_ = (id(self.func), time, RETURN)

            return RETURN



    def reconstruct(self, xi, eta, sigma, time=None, ravel=False, i=None, where='mesh-element'):
        """

        :param time:
        :param xi:
        :param eta:
        :param sigma:
        :param ravel:
        :param i:
            (1) for where == 'mesh-element' and self.ftype == "standard":
                i is None or int: the mesh element #i, if i is None, then we do it in all local mesh elements.
        :param where:
        :return:
        """
        if time is None:
            time = self.current_time
        else:
            self.current_time = time

        if where == 'mesh-element':  # input `i` means mesh element, we reconstruct it in mesh elements

            xi, eta, sigma = np.meshgrid(xi, eta, sigma, indexing='ij')
            xyz = dict()
            value = dict()

            if self.ftype == "standard":
                assert isinstance(i, int) or i is None, f"We currently only accept int or None for i"
                INDICES = self.mesh.elements.indices if i is None else [i,]
                func = self.___DO_evaluate_func_at_time___(time)
                for i in INDICES:
                    element = self.mesh.elements[i]
                    xyz_i = element.coordinate_transformation.mapping(xi, eta, sigma)
                    v00, v01, v02 = func[0][0](*xyz_i), func[0][1](*xyz_i), func[0][2](*xyz_i)
                    v10, v11, v12 = func[1][0](*xyz_i), func[1][1](*xyz_i), func[1][2](*xyz_i)
                    v20, v21, v22 = func[2][0](*xyz_i), func[2][1](*xyz_i), func[2][2](*xyz_i)
                    if ravel:
                        xyz[i] = [I.ravel('F') for I in xyz_i]
                        value[i] = ([v00.ravel('F'), v01.ravel('F'), v02.ravel('F')],
                                    [v10.ravel('F'), v11.ravel('F'), v12.ravel('F')],
                                    [v20.ravel('F'), v21.ravel('F'), v22.ravel('F')])
                    else:
                        xyz[i] = xyz_i
                        value[i] = ([v00, v01, v02],
                                    [v10, v11, v12],
                                    [v20, v21, v22])
            else:
                raise NotImplementedError(f"mesh-reconstruct not implemented for ftype: {self.ftype}")

            return xyz, value


        else:
            raise NotImplementedError(f"Can not reconstruct 3dCSCG tensor field on {where}.")




    @property
    def shape(self):
        return 3, 3

    @property
    def DO(self):
        return self._DO_

    @property
    def numerical(self):
        """The numerical property: A wrapper of all numerical methods, properties."""
        if self._numerical_ is None:
            self._numerical_ = _3dCSCG_TensorField_Numerical(self)
        return self._numerical_




class _3dCSCG_TensorField_DO(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._tf_.___DO_evaluate_func_at_time___(time=time)

    def reconstruct(self, *args, **kwargs):
        return self._tf_.reconstruct(*args, **kwargs)



class _3dCSCG_TensorField_Numerical(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    @property
    def time_derivative(self):
        """Return a _3dCSCG_TensorField instances which is the numerical time derivative of self."""
        if self._tf_.ftype == 'standard':
            F0, F1, F2 = self._tf_.func
            F00, F01, F02 = F0
            F10, F11, F12 = F1
            F20, F21, F22 = F2
            NPD4F_00 = NumericalPartialDerivative_txyz_Functions(F00)
            NPD4F_01 = NumericalPartialDerivative_txyz_Functions(F01)
            NPD4F_02 = NumericalPartialDerivative_txyz_Functions(F02)
            NPD4F_10 = NumericalPartialDerivative_txyz_Functions(F10)
            NPD4F_11 = NumericalPartialDerivative_txyz_Functions(F11)
            NPD4F_12 = NumericalPartialDerivative_txyz_Functions(F12)
            NPD4F_20 = NumericalPartialDerivative_txyz_Functions(F20)
            NPD4F_21 = NumericalPartialDerivative_txyz_Functions(F21)
            NPD4F_22 = NumericalPartialDerivative_txyz_Functions(F22)
            TDT = _3dCSCG_TensorField(self._tf_.mesh,
                                      [(NPD4F_00('t'), NPD4F_01('t'), NPD4F_02('t')),
                                       (NPD4F_10('t'), NPD4F_11('t'), NPD4F_12('t')),
                                       (NPD4F_20('t'), NPD4F_21('t'), NPD4F_22('t'))],
                                      ftype='standard',
                                      valid_time=self._tf_.valid_time,
                                      name='time-derivative-of-' + self._tf_.standard_properties.name
                                      )
            return TDT
        else:
            raise NotImplementedError(
                f"Numerical time derivative not implemented for tensor type = {self._tf_.ftype}.")

    @property
    def divergence(self):
        """Return a _3dCSCG_VectorField instances which is the divergence of self."""
        if self._tf_.ftype == 'standard':
            F0, F1, F2 = self._tf_.func
            F00, F01, F02 = F0
            F10, F11, F12 = F1
            F20, F21, F22 = F2
            NPD4F_00 = NumericalPartialDerivative_txyz_Functions(F00)
            NPD4F_01 = NumericalPartialDerivative_txyz_Functions(F01)
            NPD4F_02 = NumericalPartialDerivative_txyz_Functions(F02)
            NPD4F_10 = NumericalPartialDerivative_txyz_Functions(F10)
            NPD4F_11 = NumericalPartialDerivative_txyz_Functions(F11)
            NPD4F_12 = NumericalPartialDerivative_txyz_Functions(F12)
            NPD4F_20 = NumericalPartialDerivative_txyz_Functions(F20)
            NPD4F_21 = NumericalPartialDerivative_txyz_Functions(F21)
            NPD4F_22 = NumericalPartialDerivative_txyz_Functions(F22)
            F00_x, F01_y, F02_z = NPD4F_00('x'), NPD4F_01('y'), NPD4F_02('z')
            F10_x, F11_y, F12_z = NPD4F_10('x'), NPD4F_11('y'), NPD4F_12('z')
            F20_x, F21_y, F22_z = NPD4F_20('x'), NPD4F_21('y'), NPD4F_22('z')
            div_func_0 = ___TENSOR_DIVERGENCE_HELPER___(F00_x, F01_y, F02_z)
            div_func_1 = ___TENSOR_DIVERGENCE_HELPER___(F10_x, F11_y, F12_z)
            div_func_2 = ___TENSOR_DIVERGENCE_HELPER___(F20_x, F21_y, F22_z)
            vector_class = getattr(import_module('_3dCSCG.field.vector'), '_3dCSCG_VectorField')
            divergence_vector = vector_class(self._tf_.mesh,
                                             (div_func_0, div_func_1, div_func_2),
                                             ftype='standard',
                                             valid_time=self._tf_.valid_time,
                                             name = 'divergence-of-' + self._tf_.standard_properties.name
                                             )
            return divergence_vector
        else:
            raise NotImplementedError(f"Numerical divergence not implemented for tensor type = {self._tf_.ftype}.")



class ___TENSOR_DIVERGENCE_HELPER___(object):
    def __init__(self, fx, fy, fz):
        self._fx_ = fx
        self._fy_ = fy
        self._fz_ = fz

    def __call__(self, t, x, y, z):
        return self._fx_(t, x, y, z) + self._fy_(t, x, y, z) + self._fz_(t, x, y, z)


if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\field\tensor.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([1,1,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.cos(3*np.pi*z)

    TS = FC('tensor', ([p,0,0], [0,p,0], [0,0,p]))

    print(TS.shape)
