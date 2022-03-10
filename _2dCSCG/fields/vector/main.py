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
from _2dCSCG.fields.base.main import _2dCSCG_Continuous_FORM_BASE
from functools import partial
from screws.functions.time_plus_2d_space._0_ import _0t_
# from scipy import sparse as spspa
from _2dCSCG.fields.vector.do import _2dCSCG_VectorField_DO



class _2dCSCG_VectorField(_2dCSCG_Continuous_FORM_BASE, ndim=2):
    """The continuous vector field."""
    def __init__(self, mesh, func, ftype='standard', valid_time=None, name='vector-field'):
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_vector_field')
        self.standard_properties.name = name
        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._do_ = _2dCSCG_VectorField_DO(self)
        self._freeze_self_()

    def ___PRIVATE_set_func___(self, func, ftype='standard'):
        """
        Use this method to set up the function body and function type.

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
                else:
                    raise Exception()
            self._func_ = func
        else:
            raise Exception(f" <_2dCSCG_VectorField> do not accept funcType={ftype}")
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

        self.current_time = time
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
        return (2, )

    def reconstruct(self, xi, eta, time=None,ravel=False, i=None):
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
            func = self.___DO_evaluate_func_at_time___(time)
            for i in INDICES:
                element = self.mesh.elements[i]
                xyz_i = element.coordinate_transformation.mapping(xi, eta)
                vx_i = func[0](*xyz_i)
                vy_i = func[1](*xyz_i)

                if ravel:
                    xyz[i] = [I.ravel('F') for I in xyz_i]
                    value[i] = [vx_i.ravel('F'), vy_i.ravel('F')]
                else:
                    xyz[i] = xyz_i
                    value[i] = [vx_i, vy_i]
        else:
            raise NotImplementedError(f"reconstruct not implemented for ftype: {self.ftype}")
        return xyz, value

    @property
    def do(self):
        return self._do_





if __name__ == '__main__':
    # mpiexec -n 6 python _2dCSCG\field\vector.py
    from _2dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def p(t, x, y): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * (x+1) * (y+1)
    def q(t, x, y): return t + np.cos(2*np.pi*x) * np.cos(np.pi*y) * (x+1) * (y+1)
    VV = FC('vector', [p,q])
    VV.current_time = 1

    VV.visualize()
