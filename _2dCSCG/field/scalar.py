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
from _2dCSCG.field.main import _2dCSCG_Continuous_FORM_BASE
from functools import partial
import numpy as np
from SCREWS.functions._2d import _0t_
from SCREWS.frozen import FrozenOnly


class _2dCSCG_ScalarField(_2dCSCG_Continuous_FORM_BASE, ndim=2):
    """The continuous scalar field."""
    def __init__(self, mesh, func, ftype='standard', valid_time=None, name='scalar-field'):
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_scalar_field')
        self.standard_properties.name = name
        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._DO_ = _2dCSCG_ScalarField_DO(self)
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

    def reconstruct(self, xi, eta, time=None, ravel=False, i=None):
        """

        :param time:
        :param xi:
        :param eta:
        :param ravel:
        :param i:
        :return:
        """
        if time is None:
            time = self.current_time
        else:
            self.current_time = time

        xi, eta = np.meshgrid(xi, eta, indexing='ij')
        xyz = dict()
        value = dict()
        if self.ftype == 'standard':
            INDICES = self.mesh.elements.indices if i is None else [i,]
            func = self.___DO_evaluate_func_at_time___(time)[0]
            for i in INDICES:
                element = self.mesh.elements[i]
                xyz_i = element.coordinate_transformation.mapping(xi, eta)
                v_i = func(*xyz_i)

                if ravel:
                    xyz[i] = [I.ravel('F') for I in xyz_i]
                    value[i] = [v_i.ravel('F'),]
                else:
                    xyz[i] = xyz_i
                    value[i] = [v_i,]
        else:
            raise NotImplementedError(f"reconstruct not implemented for ftype: {self.ftype}")
        return xyz, value

    @property
    def DO(self):
        return self._DO_



class _2dCSCG_ScalarField_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._sf_.___DO_evaluate_func_at_time___(time=time)

    def reconstruct(self, *args, **kwargs):
        return self._sf_.reconstruct(*args, **kwargs)



if __name__ == '__main__':
    # mpiexec -n 6 python _2dCSCG\field\scalar.py
    from _2dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y)
    SS = FC('scalar', p)
    SS.current_time = 100

    SS.visualize()

