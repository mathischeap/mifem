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
from _3dCSCG.field.base.main import _3dCSCG_Continuous_FORM_BASE
from functools import partial
from root.config import *
from SCREWS.frozen import FrozenOnly
from SCREWS.functions._4d import CFG

from importlib import import_module
from SCREWS.numerical._4d import NumericalPartialDerivative_txyz_Functions








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

    def reconstruct(self, xi, eta, sigma, time=None, ravel=False, i=None, where=None):
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
        # we deal with default `where` input ---------------------------------------------------------------
        if where is None:
            if self.ftype == "standard":
                where = "mesh-element"
            elif self.ftype in ("boundary-wise", "trace-element-wise"):
                where = "trace-element"
            else:
                where = "mesh-element"
        else:
            pass

        # we deal with `time` input ---------------------------------------------------------------
        if time is None:
            time = self.current_time
        else:
            self.current_time = time

        # we get the current function ---------------------------------------------------------------
        func = self.___DO_evaluate_func_at_time___(time)

        # we do the reconstruction accordingly ---------------------------------------------------------------
        if where == 'mesh-element': # input `i` means mesh element, we reconstruct it in mesh elements
            xi, eta, sigma = np.meshgrid(xi, eta, sigma, indexing='ij')
            xyz = dict()
            value = dict()

            if self.ftype == "standard":
                assert isinstance(i, int) or i is None, f"We currently only accept int or None for i"
                INDICES = self.mesh.elements.indices if i is None else [i, ]

                for i in INDICES:
                    element = self.mesh.elements[i]
                    xyz_i = element.coordinate_transformation.mapping(xi, eta, sigma)
                    v_i = func[0](*xyz_i)

                    if ravel:
                        xyz[i] = [I.ravel('F') for I in xyz_i]
                        value[i] = [v_i.ravel('F'),]
                    else:
                        xyz[i] = xyz_i
                        value[i] = [v_i,]
            else:
                raise NotImplementedError(f"_3dCSCG_ScalarField mesh-reconstruct not implemented for ftype: {self.ftype}")

            return xyz, value

        elif where == 'trace-element': # input `i` means trace element, we reconstruct it on trace elements

            xyz = dict()
            value = dict()

            if self.ftype == 'boundary-wise':

                if i in (None, 'on_mesh_boundaries'):
                    RTE = self.mesh.boundaries.RANGE_trace_elements
                    INDICES = list()
                    for bn in self.func:
                        INDICES.extend(RTE[bn])

                else:
                    raise NotImplementedError(f"_3dCSCG_ScalarField of ftype 'boundary-wise'"
                                              f"trace-element-reconstruction currently doesn't accept i={i}.")

                for I in INDICES:
                    te = self.mesh.trace.elements[I]
                    assert te.IS_on_mesh_boundary, f"must be the case!"
                    xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, parse_3_1d_eps=True)

                    bn = te.on_mesh_boundary
                    assert bn in func, f"trace element #{I} is on <{bn}> which is not covered by boundary-wise func."
                    func_i = func[bn][0]

                    v_i = func_i(*xyz_i)

                    if ravel:
                        xyz[I] = [_.ravel('F') for _ in xyz_i]
                        value[I] = [v_i.ravel('F'),]
                    else:
                        xyz[I] = xyz_i
                        value[I] = [v_i,]

            elif self.ftype == 'trace-element-wise':

                if i is None: # we reconstruct on all valid local trace elements
                    INDICES = list()
                    # noinspection PyUnresolvedReferences
                    INDICES.extend(func.keys())
                elif i == 'on_mesh_boundaries': # we only reconstruct on all the valid local trace elements which are also on mesh boundaries.
                    CMB = self.covered_mesh_boundaries # will contain all mesh boundary names.
                    RTE = self.mesh.boundaries.RANGE_trace_elements
                    boundary_trace_elements = list() # local trace elements on all mesh boundaries
                    for mb in CMB:
                        boundary_trace_elements.extend(RTE[mb])
                    ___ = list()
                    # noinspection PyUnresolvedReferences
                    ___.extend(func.keys())
                    INDICES = list()
                    for I in ___:
                        if I in boundary_trace_elements:
                            INDICES.append(I)

                else:
                    raise NotImplementedError(f"_3dCSCG_ScalarField of 'trace-element-wise' ftype "
                                              f"trace-element-reconstruction currently don't accept i={i}."
                                              f"i must be one of (None, 'on_mesh_boundaries').")

                for I in INDICES: # go through all valid local trace elements

                    xyz_i, v_i = func[I](xi, eta, sigma)

                    if ravel:
                        xyz[I] = [_.ravel('F') for _ in xyz_i]
                        value[I] = [v_i[0].ravel('F') ,]
                    else:
                        xyz[I] = xyz_i
                        value[I] = v_i


            else:
                raise NotImplementedError(f"_3dCSCG_ScalarField trace-reconstruct not implemented for ftype: {self.ftype}")

            return xyz, value

        else:
            raise NotImplementedError(f"_3dCSCG_ScalarField cannot reconstruct on {where}.")



    @property
    def DO(self):
        return self._DO_

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



class ___SCALAR_NEG_HELPER_1___(object):
    def __init__(self, v):
        self._v_ = v

    def __call__(self, t, x, y, z):
        return - self._v_(t, x, y, z)

class ___SCALAR_SUB_HELPER_1___(object):
    def __init__(self, w, u):
        self._w_ = w
        self._u_ = u

    def __call__(self, t, x, y, z):
        return self._w_(t, x, y, z) - self._u_(t, x, y, z)

class ___SCALAR_ADD_HELPER_1___(object):
    def __init__(self, w, u):
        self._w_ = w
        self._u_ = u

    def __call__(self, t, x, y, z):
        return self._w_(t, x, y, z) + self._u_(t, x, y, z)





class _3dCSCG_ScalarField_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._sf_.___DO_evaluate_func_at_time___(time=time)

    def reconstruct(self, *args, **kwargs):
        return self._sf_.reconstruct(*args, **kwargs)




class _3dCSCG_ScalarField_Numerical(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    @property
    def time_derivative(self):
        """Return a _3dCSCG_ScalarField instances which is the numerical time derivative of self."""
        if self._sf_.ftype == 'standard':
            func = self._sf_.func[0]
            NPD4F = NumericalPartialDerivative_txyz_Functions(func)

            TDS = _3dCSCG_ScalarField(self._sf_.mesh, NPD4F('t'),
                              ftype='standard',
                              valid_time=self._sf_.valid_time,
                              name = 'time-derivative-of-' + self._sf_.standard_properties.name
                              )
            return TDS

        else:
            raise NotImplementedError(f"Numerical time derivative not implemented for scalar type = {self._sf_.ftype}.")

    @property
    def gradient(self):
        """Return a _3dCSCG_VectorField instances which is the numerical gradient of self."""
        if self._sf_.ftype == 'standard':
            func = self._sf_.func[0]
            NPD4F = NumericalPartialDerivative_txyz_Functions(func)
            vector_class = getattr(import_module('_3dCSCG.field.vector'), '_3dCSCG_VectorField')

            GV = vector_class(self._sf_.mesh, (NPD4F('x'), NPD4F('y'), NPD4F('z')),
                              ftype='standard',
                              valid_time=self._sf_.valid_time,
                              name = 'gradient-of-' + self._sf_.standard_properties.name
                              )
            return GV
        else:
            raise NotImplementedError(f"Numerical gradient not implemented for scalar type = {self._sf_.ftype}.")





if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\field\scalar.py
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
