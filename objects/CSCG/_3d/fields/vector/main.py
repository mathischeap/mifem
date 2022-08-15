# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
from types import FunctionType, MethodType
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_ColumnVector
from objects.CSCG._3d.fields.base import _3dCSCG_Continuous_FORM_BASE
from functools import partial
from screws.functions.time_plus_3d_space.constant import CFG

from importlib import import_module
from objects.CSCG._3d.fields.vector.numerical.main import _3dCSCG_VectorField_Numerical
from objects.CSCG._3d.fields.vector.do.main import _3dCSCG_VectorField_DO
from objects.CSCG._3d.fields.vector.components.main import _3dCSCG_VectorField_Components
from objects.CSCG._3d.fields.vector.helpers.neg import ___VECTOR_NEG_HELPER_1___
from objects.CSCG._3d.fields.vector.helpers.sub import ___VECTOR_SUB_HELPER_1___
from objects.CSCG._3d.fields.vector.helpers.add import ___VECTOR_ADD_HELPER_1___
from objects.CSCG._3d.fields.vector.helpers.flux import ___VECTOR_FLUX___
from objects.CSCG._3d.fields.vector.helpers.inner_product_with_1form import _VF_InnerWith1Form
from objects.CSCG._3d.fields.vector.helpers.inner_product_with_2form import _VF_InnerWith2Form

from objects.CSCG._3d.fields.vector.visualize.main import _3dCSCG_VectorField_Visualize


class _3dCSCG_VectorField(_3dCSCG_Continuous_FORM_BASE, ndim=3):
    """The continuous vector field."""
    def __init__(self, mesh, func, ftype=None, valid_time=None, name='vector-field'):
        if ftype is None:
            if isinstance(func, dict):
                ftype= 'boundary-wise'
            else:
                ftype = 'standard'
        else:
            pass
        super().__init__(mesh, ftype, valid_time)
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_vector_field')
        self.standard_properties.name = name
        self.___PRIVATE_set_func___(func, ftype=ftype)
        self._previous_func_id_time_ = (None, None, None)
        self._DO_ = _3dCSCG_VectorField_DO(self)
        self._visualize_ = _3dCSCG_VectorField_Visualize(self)
        self._numerical_ = None
        self._components_ = None
        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"3dCSCG_vector_field=<{self.standard_properties.name}>@{id(self)}"

    def ___PRIVATE_set_func___(self, func, ftype='standard'):
        """
        Use this method to set up the function body and function type.

        Whenever define a new funcType, edit the currentFunc for the new type.
        """
        if ftype == 'standard':
            # standard func is function or method.
            assert len(func) == 3, f"Standard vector only accepts list or tuple of shape (3,)."
            _func_checked_ = list()
            for i, fci in enumerate(func):
                if isinstance(fci, FunctionType):
                    # noinspection PyUnresolvedReferences
                    assert fci.__code__.co_argcount >= 4
                elif isinstance(fci, MethodType):
                    # noinspection PyUnresolvedReferences
                    assert fci.__code__.co_argcount >= 5
                elif isinstance(fci, (int, float)):
                    fci = CFG(fci)()
                elif callable(fci): # any other callable objects, we do not do check anymore.
                    pass
                else:
                    raise Exception(f"func[{i}]={fci} is wrong!")
                _func_checked_.append(fci)

            self._func_ = _func_checked_

        elif ftype == 'boundary-wise': # only valid (still as a vector) on mesh boundary (not domain boundary-wise)
            assert isinstance(func, dict), f" when ftype == 'boundary-wise', " \
                                           f"we must put functions in a dict whose " \
                                           f"keys are boundary names and values are" \
                                           f"the functions."
            for bn in func:
                assert bn in self.mesh.boundaries.names, \
                    f"func key: [{bn}] is not a valid boundary name " \
                    f"({self.mesh.boundaries.names})"

                func_bn = func[bn]
                assert len(func_bn) == 3, \
                    f"3d vector should be of shape (3,), now it is {np.shape(func_bn)}."

                _func_bn_ck_ = list()
                for fci in func_bn:
                    # standard func is function or method.
                    if isinstance(fci, FunctionType):
                        assert fci.__code__.co_argcount >= 4
                    elif isinstance(fci, MethodType):
                        # noinspection PyUnresolvedReferences
                        assert fci.__code__.co_argcount >= 5
                    elif isinstance(fci, (int, float)):
                        fci = CFG(fci)()
                    else:
                        raise Exception()
                    _func_bn_ck_.append(fci)

                func[bn] = _func_bn_ck_

            self._func_ = func

        elif ftype == 'trace-element-wise':
            # we have received a dict whose keys are local trace elements, values are callable that returns, xyz and a vector.
            assert isinstance(func, dict), f"func for trace-element-wise vector must a dict."
            for i in func: # valid local trace elements
                assert i in self.mesh.trace.elements, f"trace element #{i} is not in this core (#{rAnk})."
                # NOTE that we do not put the vector in a list or tuple, it should take (t, xi, eta, sigma) and then return xyz and the vector.
                assert callable(func[i]), f"func[{i}] is not callable."
            self._func_ = func
        else:
            raise Exception(f" <_3dCSCG_VectorField> do not accept funcType={ftype}")

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
                RETURN = partial(self.func[0], time), \
                         partial(self.func[1], time), \
                         partial(self.func[2], time)

            elif self.ftype  == 'boundary-wise':

                RETURN = dict()
                for bn in self.func:
                    RETURN[bn] = [partial(self.func[bn][0], time),
                                  partial(self.func[bn][1], time),
                                  partial(self.func[bn][2], time)]

            elif self.ftype  == 'trace-element-wise':
                RETURN = dict()
                for i in self.func: # go through all valid trace elements
                    vi = self.func[i]
                    RETURN[i] = partial(vi, time) # We can see that for each trace-element, it is a single function

            else:
                raise Exception(f" do not understand funcType={self.ftype}")


            self._previous_func_id_time_ = (id(self.func), time, RETURN)

            return RETURN

    @property
    def shape(self):
        return 3, # do not remove comma.

    def reconstruct(self, *args, **kwargs):
        return self.do.reconstruct(*args, **kwargs)


    def ___PRIVATE_do_inner_product_with_space_of___(self, other, quad_degree=None):
        """
        do :math:`(\\cdot, \\cdot)` with the basis functions of given standard form.

        The time will be the current time.

        :param other:
        :param quad_degree:
        """
        if self.ftype == 'standard':
            assert self.mesh == other.mesh
            if other.__class__.__name__ == '_2Form':
                # when other is a 2-form, we consider self as a (star of) continuous 1-form.
                DG = _VF_InnerWith2Form(self, other, quad_degree)
            elif other.__class__.__name__ == '_1Form':
                # when other is a 1-form, we consider self as a (star of) continuous 2-form.
                DG = _VF_InnerWith1Form(self, other, quad_degree)
            else:
                raise Exception(f"standard _3dCSCG_VectorField can only inner-product with 1- or 2-standard form.")
            return EWC_ColumnVector(self.mesh.elements, DG)
        else:
            raise NotImplementedError(f"_3dCSCG_VectorField of ftype='{self.ftype}' "
                                      f"cannot inner product with {other.__class__}")

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
            self._numerical_ = _3dCSCG_VectorField_Numerical(self)
        return self._numerical_

    @property
    def components(self):
        """A wrapper of all components of this vector"""
        if self._components_ is None:
            self._components_ = _3dCSCG_VectorField_Components(self)
        return self._components_

    @property
    def flux(self):
        """Return a _3dCSCG_ScalarField representing the flux scalar on all valid trace elements.

        Let the self vector is u, then we return a scalar (u \dot n) where n is the positive unit norm vector.

        When the self vector is of ftype
            - 'standard': we will make a ('trace-element-wise') scalar valid on all local trace elements.

        """
        if self.ftype == 'standard':
            # we have a standard vector, we will make a flux scalar valid on all (locally in each core) trace elements.

            safe_copy = _3dCSCG_VectorField(self.mesh,
                                            self.func,
                                            ftype=self.ftype,
                                            valid_time=self.valid_time,
                                            name=self.standard_properties.name
                                            ) # we made a safe copy.
            # this is very important as it decoupled the norm component and the vector. MUST do THIS!

            trace_element_wise_func = dict()
            for i in safe_copy.mesh.trace.elements: # the local trace element #i on mesh boundaries
                trace_element_wise_func[i] = ___VECTOR_FLUX___(safe_copy, i)

            base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-3]) + '.'

            scalar_class = getattr(import_module(base_path + 'scalar.main'), '_3dCSCG_ScalarField')
            return scalar_class(safe_copy.mesh, trace_element_wise_func,
                                       ftype='trace-element-wise',
                                       valid_time=safe_copy.valid_time,
                                       name='flux-scalar-of-' + safe_copy.standard_properties.name
                                       )

        else:
            raise NotImplementedError(f"`flux` scalar of a vector of "
                                      f"type = {self.ftype} is not implemented")

    def __neg__(self):
        """-self."""
        if self.ftype == 'standard':
            w0, w1, w2 = self.func

            x0 = ___VECTOR_NEG_HELPER_1___(w0)
            x1 = ___VECTOR_NEG_HELPER_1___(w1)
            x2 = ___VECTOR_NEG_HELPER_1___(w2)

            neg_vector = _3dCSCG_VectorField(self.mesh,
                 [x0, x1, x2],
                 ftype='standard',
                 valid_time=self.valid_time,
                 name = '-' + self.standard_properties.name
                                            )
            return neg_vector

        else:
            raise Exception(f"cannot do neg for {self.ftype} _3dCSCG_VectorField.")

    def __sub__(self, other):
        """self - other"""
        if other.__class__.__name__ == '_3dCSCG_VectorField':

            if self.ftype == 'standard' and other.ftype == 'standard':

                w0, w1, w2 = self.func
                u0, u1, u2 = other.func

                x0 = ___VECTOR_SUB_HELPER_1___(w0, u0)
                x1 = ___VECTOR_SUB_HELPER_1___(w1, u1)
                x2 = ___VECTOR_SUB_HELPER_1___(w2, u2)

                sub_vector = _3dCSCG_VectorField(self.mesh,
                     [x0, x1, x2],
                     ftype='standard',
                     valid_time=self.valid_time,
                     name = self.standard_properties.name + '-' + other.standard_properties.name
                                                )
                return sub_vector

            else:
                raise Exception(f"cannot do {self.ftype} _3dCSCG_VectorField - {other.ftype} _3dCSCG_VectorField")
        else:
            raise Exception(f"cannot do _3dCSCG_VectorField - {other.__class__}")

    def __add__(self, other):
        """self + other"""
        if other.__class__.__name__ == '_3dCSCG_VectorField':

            if self.ftype == 'standard' and other.ftype == 'standard':

                w0, w1, w2 = self.func
                u0, u1, u2 = other.func

                x0 = ___VECTOR_ADD_HELPER_1___(w0, u0)
                x1 = ___VECTOR_ADD_HELPER_1___(w1, u1)
                x2 = ___VECTOR_ADD_HELPER_1___(w2, u2)

                add_vector = _3dCSCG_VectorField(self.mesh,
                     [x0, x1, x2],
                     ftype='standard',
                     valid_time=self.valid_time,
                     name = self.standard_properties.name + '+' + other.standard_properties.name
                                                )
                return add_vector

            else:
                raise Exception(f"cannot do {self.ftype} _3dCSCG_VectorField - {other.ftype} _3dCSCG_VectorField")
        else:
            raise Exception(f"cannot do _3dCSCG_VectorField + {other.__class__}")






if __name__ == '__main__':
    # mpiexec -n 6 python objects/CSCG/_3d/fields/vector/main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([1,1,1], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    def velocity_x(t, x, y, z): return t + np.cos(2*np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z)
    def velocity_y(t, x, y, z): return t + np.sin(np.pi*x) * np.cos(2*np.pi*y) * np.sin(np.pi*z)
    def velocity_z(t, x, y, z): return t + np.sin(np.pi*x) * np.sin(np.pi*y) * np.cos(2*np.pi*z)

    SV = FC('vector', [velocity_x, velocity_y, velocity_z])


    norm = SV.components.norm
    norm.current_time=0
    norm.visualize()

    para = SV.components.T_para
    para.current_time=0
    para.visualize()

    perp = SV.components.T_perp
    perp.current_time=0
    perp.visualize()