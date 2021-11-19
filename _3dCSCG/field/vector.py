# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
from types import FunctionType, MethodType
from SCREWS.frozen import FrozenOnly
from TOOLS.linear_algebra.elementwise_cache import EWC_ColumnVector
from _3dCSCG.field.main import _3dCSCG_Continuous_FORM_BASE
from functools import partial
from scipy import sparse as spspa
from SCREWS.functions._4d import CFG

from importlib import import_module
from SCREWS.numerical._4d import NumericalPartialDerivative_txyz_Functions

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
        self._numerical_ = None
        self._components_ = None
        self._freeze_self_()

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
                RETURN = partial(self.func[0], time), partial(self.func[1], time), partial(self.func[2], time)

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
        return (3,)

    def reconstruct(self, xi, eta, sigma,
        time=None,
        ravel=False,
        i=None, where='mesh-element',
        structured=True):
        """

        :param xi: if `structured`, `xi` must be 1d array.
        :param eta: if `structured`, `eta` must be 1d array.
        :param sigma: if `structured`, `sigma` must be 1d array.
        :param time:
        :param ravel:
        :param i:
            (1) for where == 'mesh-element' and self.ftype == "standard":
                i is None or int: the mesh element #i, if i is None, then we do it in all local mesh elements.
        :param where:
        :param structured: When `structured`, we must have xi, eta and sigma be 1d, so we do meshgrid from
            them. This is to keep it consistent with the reconstruct of forms.

            When not `structured`, then it is free. we reconstruct on (xi, eta, sigma), no matter the dimensions.

        :return:
        """
        if time is None:
            time = self.current_time
        else:
            self.current_time = time


        if structured:
            pass # we stay at this method
        else:
            # we go to the unstructured reconstruction method and return here!
            return self.___PRIVATE_unstructured_reconstruction___(
                xi, eta, sigma, time, ravel, i, where, structured)


        func = self.___DO_evaluate_func_at_time___(time)

        if where == 'mesh-element':  # input `i` means mesh element, we reconstruct it in mesh elements

            xi, eta, sigma = np.meshgrid(xi, eta, sigma, indexing='ij')
            xyz = dict()
            value = dict()


            if self.ftype == "standard":
                if isinstance(i, int):
                    INDICES = [i,]
                elif i is None:
                    INDICES = self.mesh.elements.indices

                else:
                    raise NotImplementedError(f"_3dCSCG_VectorField of 'standard' ftype"
                                              f" mesh-element-reconstruction currently doesn't accept i={i}.")

                for i in INDICES:
                    element = self.mesh.elements[i]
                    xyz_i = element.coordinate_transformation.mapping(xi, eta, sigma)
                    vx_i = func[0](*xyz_i)
                    vy_i = func[1](*xyz_i)
                    vz_i = func[2](*xyz_i)

                    if ravel:
                        xyz[i] = [I.ravel('F') for I in xyz_i]
                        value[i] = [vx_i.ravel('F'), vy_i.ravel('F'), vz_i.ravel('F')]
                    else:
                        xyz[i] = xyz_i
                        value[i] = [vx_i, vy_i, vz_i]
            else:
                raise NotImplementedError(f"mesh-reconstruct not implemented for ftype: {self.ftype}")

            return xyz, value

        elif where == 'trace-element': # input `i` means trace element, we reconstruct it in mesh elements

            xyz = dict()
            value = dict()

            if self.ftype == 'standard':
                if isinstance(i, int):
                    INDICES = [i,]
                else:
                    raise NotImplementedError(f"_3dCSCG_VectorField of 'standard' ftype"
                                              f" trace-element-reconstruction currently don't accept i={i}.")

                for I in INDICES:
                    te = self.mesh.trace.elements[I]
                    xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, parse_3_1d_eps=True)

                    vx_i = func[0](*xyz_i)
                    vy_i = func[1](*xyz_i)
                    vz_i = func[2](*xyz_i)

                    if ravel:
                        xyz[I] = [_.ravel('F') for _ in xyz_i]
                        value[I] = [vx_i.ravel('F'), vy_i.ravel('F'), vz_i.ravel('F')]
                    else:
                        xyz[I] = xyz_i
                        value[I] = [vx_i, vy_i, vz_i,]

            elif self.ftype == 'boundary-wise':


                if i in (None, 'on_mesh_boundaries'):
                    INDICES = list()
                    RTE = self.mesh.boundaries.RANGE_trace_elements
                    for bn in self.func: # this may not contain all mesh boundaries, only valid ones.
                        INDICES.extend(RTE[bn])

                else:
                    raise NotImplementedError(f"_3dCSCG_VectorField of 'boundary-wise' ftype"
                                              f" trace-element-reconstruction currently don't accept i={i}.")

                for I in INDICES:

                    te = self.mesh.trace.elements[I]
                    assert te.IS_on_mesh_boundary, f"must be the case because ftype == 'boundary-wise!"
                    xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, parse_3_1d_eps=True)

                    bn = te.on_mesh_boundary
                    assert bn in func, f"trace element #{I} is on <{bn}> which is not covered by boundary-wise func."
                    func_i = func[bn]

                    vx_i = func_i[0](*xyz_i)
                    vy_i = func_i[1](*xyz_i)
                    vz_i = func_i[2](*xyz_i)

                    if ravel:
                        xyz[I] = [_.ravel('F') for _ in xyz_i]
                        value[I] = [vx_i.ravel('F'), vy_i.ravel('F'), vz_i.ravel('F')]
                    else:
                        xyz[I] = xyz_i
                        value[I] = [vx_i, vy_i, vz_i,]

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
                    raise NotImplementedError(f"_3dCSCG_VectorField of 'trace-element-wise' ftype "
                                              f"trace-element-reconstruction currently don't accept i={i}."
                                              f"i must be one of (None, 'on_mesh_boundaries').")

                for I in INDICES: # go through all valid local trace elements

                    xyz_i, v_i = func[I](xi, eta, sigma)

                    if ravel:
                        xyz[I] = [_.ravel('F') for _ in xyz_i]
                        value[I] = [_.ravel('F') for _ in v_i]
                    else:
                        xyz[I] = xyz_i
                        value[I] = v_i

            else:
                raise NotImplementedError(f"_3dCSCG_VectorField trace-element-wise-reconstruction "
                                          f"not implemented for ftype: {self.ftype}")

            return xyz, value

        else:
            raise NotImplementedError(f"_3dCSCG_VectorField cannot reconstruct on {where}.")

    def ___PRIVATE_unstructured_reconstruction___(self, xi, eta, sigma, time, ravel, i, where, structured):
        """
        Note that even (xi, eta, sigma) are unstructured (we will not meshgrid them), their dimension can be
        arbitrary but their shapes must be same. So we can still do ravel to the results, but the raveled
        results may be scatter randomly.

        :param xi:
        :param eta:
        :param sigma:
        :param time:
        :param ravel:
        :param i:
            (1) for where == 'mesh-element' and self.ftype == "standard":
                i is None or int: the mesh element #i, if i is None, then we do it in all local mesh elements.
        :param where:
        :param structured: it must be False

        :return:
        """
        assert not structured, f"Must be unstructured reconstruction!"

        if where == 'trace-element':  # input `i` means trace element, we reconstruct it in mesh elements

            xyz = dict()
            value = dict()
            func = self.___DO_evaluate_func_at_time___(time)

            if self.ftype == 'standard':
                RTE = self.mesh.boundaries.RANGE_trace_elements

                if i is None:
                    INDICES = list()  # will reconstruct for all local trace elements on the mesh boundaries.
                    for bn in RTE:
                        INDICES.extend(RTE[bn])
                elif isinstance(i, int):# when i is int, we reconstruct on this particular trace element #i.
                    # so this may be running locally. Cause each trace element is a local instance.
                    INDICES = [i,]
                else:
                    raise NotImplementedError(
                        f"currently only accept i=None, so reconstruct for all valid trace-elements.")

                for i in INDICES:
                    te = self.mesh.trace.elements[i]
                    assert te.IS_on_mesh_boundary, f"must be the case!"
                    xyz_i = te.coordinate_transformation.mapping(xi, eta, sigma, picking=True)

                    vx_i = func[0](*xyz_i)
                    vy_i = func[1](*xyz_i)
                    vz_i = func[2](*xyz_i)

                    if ravel:
                        xyz[i] = [I.ravel('F') for I in xyz_i]
                        value[i] = [vx_i.ravel('F'), vy_i.ravel('F'), vz_i.ravel('F')]
                    else:
                        xyz[i] = xyz_i
                        value[i] = [vx_i, vy_i, vz_i, ]

            else:
                raise NotImplementedError(f"trace-unstructured-reconstruct not implemented for ftype: {self.ftype}")

            return xyz, value

        else:
            raise NotImplementedError(f"Can not reconstruct (unstructured) on {where}.")






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
    def DO(self):
        return self._DO_

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
            # this is very important as it decoupled the norm component and the vector. MUST DO THIS!

            trace_element_wise_func = dict()
            for i in safe_copy.mesh.trace.elements: # the local trace element #i on mesh boundaries
                trace_element_wise_func[i] = ___VECTOR_FLUX___(safe_copy, i)
            scalar_class = getattr(import_module('_3dCSCG.field.scalar'), '_3dCSCG_ScalarField')
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



class ___VECTOR_NEG_HELPER_1___(object):
    def __init__(self, v):
        self._v_ = v

    def __call__(self, t, x, y, z):
        return - self._v_(t, x, y, z)

class ___VECTOR_SUB_HELPER_1___(object):
    def __init__(self, w, u):
        self._w_ = w
        self._u_ = u

    def __call__(self, t, x, y, z):
        return self._w_(t, x, y, z) - self._u_(t, x, y, z)

class ___VECTOR_ADD_HELPER_1___(object):
    def __init__(self, w, u):
        self._w_ = w
        self._u_ = u

    def __call__(self, t, x, y, z):
        return self._w_(t, x, y, z) + self._u_(t, x, y, z)



class ___VECTOR_FLUX___(object):
    """Here we will wrap the reconstruction of standard vector such that it works like a function. Then we
    can use it to build flux scalar."""
    def __init__(self, vf, i):
        """

        :param vf: the vector
        :param i: for trace element #i
        """
        self._vf_ = vf
        self._i_ = i
        self._te_ = vf.mesh.trace.elements[i]

    def __call__(self, t, xi, et, sg):
        """This actually includes a reconstruction. So we will call from xi, et, sg.

        :param t: at time t
        :param xi: must be a 1d array.
        :param et: must be a 1d array.
        :param sg: must be a 1d array.
        :return:
        """
        vf = self._vf_
        i = self._i_

        if vf.ftype == 'standard':

            xyz, w = vf.reconstruct(xi, et, sg,
                                    time=t,
                                    i=i,
                                    ravel=False,
                                    where='trace-element',
                                    structured=True)
            xyz = xyz[i]
            w = w[i]
            n = self._te_.coordinate_transformation.unit_normal_vector(xi, et, sg, parse_3_1d_eps=True)
            w_dot_n = w[0]*n[0] + w[1]*n[1] + w[2]*n[2]
            return xyz, (w_dot_n,) # xyz and the norm scalar.

        else:
            raise Exception(f"We cannot reconstruct from txyz for {self._vf_.ftype} vector.")





class _3dCSCG_VectorField_DO(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._vf_.___DO_evaluate_func_at_time___(time=time)

    def inner_product_with_space_of(self, other, quad_degree=None):
        return self._vf_.___PRIVATE_do_inner_product_with_space_of___(other, quad_degree=quad_degree)

    def reconstruct(self, *args, **kwargs):
        return self._vf_.reconstruct(*args, **kwargs)

    def cross_product(self, vector):
        """We compute self cross_product another vector: self X vector"""

        if vector.__class__.__name__ == '_3dCSCG_VectorField':

            if self._vf_.ftype == 'standard':
                if vector.ftype == 'standard':
                    w0, w1, w2 = self._vf_.func
                    u0, u1, u2 = vector.func
                    CP0 = ___VF_CROSS_PRODUCT_HELPER_1___(w1, u2, w2, u1)
                    CP1 = ___VF_CROSS_PRODUCT_HELPER_1___(w2, u0, w0, u2)
                    CP2 = ___VF_CROSS_PRODUCT_HELPER_1___(w0, u1, w1, u0)
                    cp_vector = _3dCSCG_VectorField(self._vf_.mesh,
                                                    [CP0, CP1, CP2],
                                                    ftype='standard',
                                                    valid_time=self._vf_.valid_time,
                                                    name = self._vf_.standard_properties.name
                                                         + '--cross-X-product--'
                                                         + vector.standard_properties.name
                                                    )
                    return cp_vector
                else:
                    raise NotImplementedError(
                        f"a standard _3dCSCG_VectorField cannot do cross product with a _3dCSCG_VectorField of ftype {vector.ftype}.")
            else:
                raise NotImplementedError(f"a _3dCSCG_VectorField of ftype {self._vf_.ftype} cannot do cross product.")
        else:
            raise NotImplementedError(f"a _3dCSCG_VectorField can not cross product a {vector}.")




class _VF_InnerWith2Form(FrozenOnly):
    def __init__(self, vf, _2f, quad_degree):
        if quad_degree is None: quad_degree = [_2f.dqp[i]+1 for i in range(3)]
        quad_nodes, _, quad_weights = _2f.space.DO_evaluate_quadrature(quad_degree)
        _, bf2 = _2f.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        self._g0_, self._g1_, self._g2_ = bf2
        self._JM_ = _2f.mesh.elements.coordinate_transformation.QUAD_1d.Jacobian_matrix(quad_degree, 'Gauss')
        self._mapping_ = _2f.mesh.elements.coordinate_transformation.QUAD_1d.mapping(quad_degree, 'Gauss')
        self._vf_ = vf
        self._qw_ = quad_weights
        self._mesh_ = _2f.mesh
        self._freeze_self_()

    def __call__(self, i):
        """
        :param i: # element.
        :return:
        """
        mark = self._mesh_.elements[i].type_wrt_metric.mark
        xyz = self._mapping_[i]
        _f0_, _f1_, _f2_ = self._vf_.DO.evaluate_func_at_time()
        f0, f1, f2 = _f0_(*xyz), _f1_(*xyz), _f2_(*xyz)
        g0, g1, g2 = self._g0_, self._g1_, self._g2_
        JM = self._JM_[i]
        if isinstance(mark, str) and mark[:4] == 'Orth':
            v0 = np.einsum('w, iw -> i', f0 * JM[0][0] * self._qw_, g0, optimize='greedy')
            v1 = np.einsum('w, iw -> i', f1 * JM[1][1] * self._qw_, g1, optimize='greedy')
            v2 = np.einsum('w, iw -> i', f2 * JM[2][2] * self._qw_, g2, optimize='greedy')
        else:
            v0 = np.einsum('w, iw -> i', (f0*JM[0][0] + f1*JM[1][0] + f2*JM[2][0])*self._qw_, g0, optimize='greedy')
            v1 = np.einsum('w, iw -> i', (f0*JM[0][1] + f1*JM[1][1] + f2*JM[2][1])*self._qw_, g1, optimize='greedy')
            v2 = np.einsum('w, iw -> i', (f0*JM[0][2] + f1*JM[1][2] + f2*JM[2][2])*self._qw_, g2, optimize='greedy')
        RETURN = spspa.csr_matrix(np.concatenate([v0, v1, v2])).T


        return RETURN


class _VF_InnerWith1Form(FrozenOnly):
    def __init__(self, vf, _1f, quad_degree):
        if quad_degree is None: quad_degree = [_1f.dqp[i]+1 for i in range(3)]
        quad_nodes, _, quad_weights = _1f.space.DO_evaluate_quadrature(quad_degree)
        _, bf1 = _1f.DO.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        self._g0_, self._g1_, self._g2_ = bf1
        self._JM_ = _1f.mesh.elements.coordinate_transformation.QUAD_1d.Jacobian_matrix(quad_degree, 'Gauss')
        self._mapping_ = _1f.mesh.elements.coordinate_transformation.QUAD_1d.mapping(quad_degree, 'Gauss')
        self._qw_ = quad_weights
        self._mesh_ = _1f.mesh
        self._vf_ = vf
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self._J_cache_ = dict()

    def ___PRIVATE_J___(self, i, mark):
        """"""
        if mark in self._J_cache_:
            return self._J_cache_[mark]
        else:
            JM = self._JM_[i]
            if isinstance(mark, str) and mark[:4] == 'Orth':
                J00 = JM[1][1] * JM[2][2]
                J01 = None
                J02 = None
                J10 = None
                J11 = JM[2][2] * JM[0][0]
                J12 = None
                J20 = None
                J21 = None
                J22 = JM[0][0] * JM[1][1]
            else:
                J00 = JM[1][1]*JM[2][2] - JM[1][2]*JM[2][1]
                J01 = JM[2][1]*JM[0][2] - JM[2][2]*JM[0][1]
                J02 = JM[0][1]*JM[1][2] - JM[0][2]*JM[1][1]
                J10 = JM[1][2]*JM[2][0] - JM[1][0]*JM[2][2]
                J11 = JM[2][2]*JM[0][0] - JM[2][0]*JM[0][2]
                J12 = JM[0][2]*JM[1][0] - JM[0][0]*JM[1][2]
                J20 = JM[1][0]*JM[2][1] - JM[1][1]*JM[2][0]
                J21 = JM[2][0]*JM[0][1] - JM[2][1]*JM[0][0]
                J22 = JM[0][0]*JM[1][1] - JM[0][1]*JM[1][0]
            J = (J00, J01, J02, J10, J11, J12, J20, J21, J22)
            self._J_cache_[mark] = J
            return J

    def __call__(self, i):
        """
        :param i: # element.
        :return:
        """
        mark = self._mesh_.elements[i].type_wrt_metric.mark
        xyz = self._mapping_[i]
        _f0_, _f1_, _f2_ = self._vf_.DO.evaluate_func_at_time()
        f0, f1, f2 = _f0_(*xyz), _f1_(*xyz), _f2_(*xyz)
        g0, g1, g2 = self._g0_, self._g1_, self._g2_
        J00, J01, J02, J10, J11, J12, J20, J21, J22 = self.___PRIVATE_J___(i, mark)
        if isinstance(mark, str) and mark[:4] == 'Orth':
            v0 = np.einsum('w, iw -> i', f0 * J00 * self._qw_, g0, optimize='greedy')
            v1 = np.einsum('w, iw -> i', f1 * J11 * self._qw_, g1, optimize='greedy')
            v2 = np.einsum('w, iw -> i', f2 * J22 * self._qw_, g2, optimize='greedy')
        else:
            v0 = np.einsum('w, iw -> i', (f0*J00 + f1*J01 + f2*J02) * self._qw_, g0, optimize='greedy')
            v1 = np.einsum('w, iw -> i', (f0*J10 + f1*J11 + f2*J12) * self._qw_, g1, optimize='greedy')
            v2 = np.einsum('w, iw -> i', (f0*J20 + f1*J21 + f2*J22) * self._qw_, g2, optimize='greedy')
        RETURN = spspa.csr_matrix(np.concatenate([v0, v1, v2])).T

        return RETURN



class ___VF_CROSS_PRODUCT_HELPER_1___(object):
    def __init__(self, f0, f1, f2, f3):
        self._f0_ = f0
        self._f1_ = f1
        self._f2_ = f2
        self._f3_ = f3

    def __call__(self, t, x, y, z):
        return self._f0_(t, x, y, z) * self._f1_(t, x, y, z) - self._f2_(t, x, y, z) * self._f3_(t, x, y, z)


class _3dCSCG_VectorField_Numerical(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    @property
    def time_derivative(self):
        """Return a _3dCSCG_VectorField instances which is the numerical time derivative of self."""
        if self._vf_.ftype == 'standard':
            func_x, func_y, func_z = self._vf_.func
            NPD4F_x = NumericalPartialDerivative_txyz_Functions(func_x)
            NPD4F_y = NumericalPartialDerivative_txyz_Functions(func_y)
            NPD4F_z = NumericalPartialDerivative_txyz_Functions(func_z)
            TDV = _3dCSCG_VectorField(self._vf_.mesh,
                                      (NPD4F_x('t'), NPD4F_y('t'), NPD4F_z('t')),
                                      ftype='standard',
                                      valid_time=self._vf_.valid_time,
                                      name='time-derivative-of-' + self._vf_.standard_properties.name
                                      )
            return TDV
        else:
            raise NotImplementedError(
                f"Numerical time derivative not implemented for vector type = {self._vf_.ftype}.")

    @property
    def gradient(self):
        """Return a _3dCSCG_TensorField instances which is the numerical gradient of self."""
        if self._vf_.ftype == 'standard':
            func_x, func_y, func_z = self._vf_.func
            NPD4F_x = NumericalPartialDerivative_txyz_Functions(func_x)
            NPD4F_y = NumericalPartialDerivative_txyz_Functions(func_y)
            NPD4F_z = NumericalPartialDerivative_txyz_Functions(func_z)
            T00, T01, T02 = NPD4F_x('x'), NPD4F_x('y'), NPD4F_x('z')
            T10, T11, T12 = NPD4F_y('x'), NPD4F_y('y'), NPD4F_y('z')
            T20, T21, T22 = NPD4F_z('x'), NPD4F_z('y'), NPD4F_z('z')

            tensor_class = getattr(import_module('_3dCSCG.field.tensor'), '_3dCSCG_TensorField')
            gradient_tensor = tensor_class(self._vf_.mesh,
                                     [(T00, T01, T02),
                                      (T10, T11, T12),
                                      (T20, T21, T22),],
                                     ftype='standard',
                                     valid_time=self._vf_.valid_time,
                                     name = 'gradient-of-' + self._vf_.standard_properties.name
                                     )
            return gradient_tensor
        else:
            raise NotImplementedError(f"Numerical gradient not implemented for vector type = {self._vf_.ftype}.")

    @property
    def curl(self):
        """Return a _3dCSCG_TensorField instances which is the numerical curl of self."""
        if self._vf_.ftype == 'standard':
            func_x, func_y, func_z = self._vf_.func
            NPD4F_x = NumericalPartialDerivative_txyz_Functions(func_x)
            NPD4F_y = NumericalPartialDerivative_txyz_Functions(func_y)
            NPD4F_z = NumericalPartialDerivative_txyz_Functions(func_z)
            u_y, u_z = NPD4F_x('y'), NPD4F_x('z')
            v_x, v_z = NPD4F_y('x'), NPD4F_y('z')
            w_x, w_y = NPD4F_z('x'), NPD4F_z('y')

            curl_vector_0 = ___VECTOR_CURL_HELPER___(w_y, v_z)
            curl_vector_1 = ___VECTOR_CURL_HELPER___(u_z, w_x)
            curl_vector_2 = ___VECTOR_CURL_HELPER___(v_x, u_y)

            curl_vector = _3dCSCG_VectorField(self._vf_.mesh,
                                     [curl_vector_0, curl_vector_1, curl_vector_2],
                                     ftype='standard',
                                     valid_time=self._vf_.valid_time,
                                     name = 'curl-of-' + self._vf_.standard_properties.name
                                     )
            return curl_vector
        else:
            raise NotImplementedError(f"Numerical curl not implemented for vector type = {self._vf_.ftype}.")

    @property
    def divergence(self):
        """Return a _3dCSCG_ScalarField instances which is the numerical divergence of self."""
        if self._vf_.ftype == 'standard':
            func_x, func_y, func_z = self._vf_.func
            NPD4F_x = NumericalPartialDerivative_txyz_Functions(func_x)
            NPD4F_y = NumericalPartialDerivative_txyz_Functions(func_y)
            NPD4F_z = NumericalPartialDerivative_txyz_Functions(func_z)
            u_x = NPD4F_x('x')
            v_y = NPD4F_y('y')
            w_z = NPD4F_z('z')
            div_func = ___VECTOR_DIVERGENCE_HELPER___(u_x, v_y, w_z)
            scalar_class = getattr(import_module('_3dCSCG.field.scalar'), '_3dCSCG_ScalarField')
            divergence_scalar = scalar_class(self._vf_.mesh,
                                             div_func,
                                             ftype='standard',
                                             valid_time=self._vf_.valid_time,
                                             name = 'divergence-of-' + self._vf_.standard_properties.name
                                             )
            return divergence_scalar
        else:
            raise NotImplementedError(f"Numerical divergence not implemented for vector type = {self._vf_.ftype}.")




class ___VECTOR_CURL_HELPER___(object):
    def __init__(self, f0, f1):
        self._f0_ = f0
        self._f1_ = f1

    def __call__(self, t, x, y, z):
        return self._f0_(t, x, y, z) - self._f1_(t, x, y, z)


class ___VECTOR_DIVERGENCE_HELPER___(object):
    def __init__(self, f0, f1, f2):
        self._f0_ = f0
        self._f1_ = f1
        self._f2_ = f2

    def __call__(self, t, x, y, z):
        return self._f0_(t, x, y, z) + self._f1_(t, x, y, z) + self._f2_(t, x, y, z)






class _3dCSCG_VectorField_Components(FrozenOnly):
    """A wrapper of all components of the vector."""
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    @property
    def norm(self):
        """The norm component of the vector on the trace elements. Note that we consider the positive norm
        direction of each trace element instead of the outward norm direction of a mesh element.
        Note that the norm component (a vector) will be of ftype 'trace-element-wise'.

        Let the vector is u, then we return a vector (u \dot n)n where n is the positive unit norm vector.

        When the vector is of ftype
            - 'standard': we will make another ('trace-element-wise') vector valid on all trace elements.

        """
        if self._vf_.ftype == 'standard':
            # we have a standard vector, we will make a norm vector valid on all (locally in each core) trace elements.

            safe_copy = _3dCSCG_VectorField(self._vf_.mesh,
                                            self._vf_.func,
                                            ftype=self._vf_.ftype,
                                            valid_time=self._vf_.valid_time,
                                            name=self._vf_.standard_properties.name
                                            )
            # this is very important as it decoupled the norm component and the vector. MUST DO THIS!

            trace_element_wise_func = dict()
            for i in safe_copy.mesh.trace.elements: # the local trace element #i on mesh boundaries
                trace_element_wise_func[i] = ___VECTOR_NORM_COMPONENT___(safe_copy, i)

            return _3dCSCG_VectorField(safe_copy.mesh, trace_element_wise_func,
                                       ftype='trace-element-wise',
                                       valid_time=safe_copy.valid_time,
                                       name='norm-component-of-' + safe_copy.standard_properties.name
                                       )

        else:
            raise NotImplementedError(f"`norm` component of a vector of "
                                      f"type = {self._vf_.ftype} is not implemented")

    @property
    def T_para(self):
        """The trace_parallel component of the vector.

        Note that the component (a vector) will be of ftype 'trace-element-wise'.

        When the vector is of ftype
            - 'standard': we will make another ('trace-element-wise') vector valid on all trace elements.
        """
        if self._vf_.ftype == 'standard':
            # we have a standard vector, we will make a T_para vector valid on all (locally in each core) trace elements.

            safe_copy = _3dCSCG_VectorField(self._vf_.mesh,
                                            self._vf_.func,
                                            ftype=self._vf_.ftype,
                                            valid_time=self._vf_.valid_time,
                                            name=self._vf_.standard_properties.name
                                            )
            # this is very important as it decoupled the norm component and the vector. MUST DO THIS!

            trace_element_wise_func = dict()
            for i in safe_copy.mesh.trace.elements: # the local trace element #i on mesh boundaries
                trace_element_wise_func[i] = ___VECTOR_T_PARA_COMPONENT___(safe_copy, i)

            return _3dCSCG_VectorField(safe_copy.mesh, trace_element_wise_func,
                                       ftype='trace-element-wise',
                                       valid_time=safe_copy.valid_time,
                                       name='T-para-component-of-' + safe_copy.standard_properties.name
                                       )

        else:
            raise NotImplementedError(f"`T_para` component of a vector of "
                                      f"type = {self._vf_.ftype} is not implemented")

    @property
    def T_perp(self):
        """The trace_perpendicular component of the vector.

        Note that the component (a vector) will be of ftype 'trace-element-wise'.

        When the vector is of ftype
            - 'standard': we will make another ('trace-element-wise') vector valid on all trace elements.
        """
        if self._vf_.ftype == 'standard':
            # we have a standard vector, we will make a T_perp vector valid on all (locally in each core) trace elements.

            safe_copy = _3dCSCG_VectorField(self._vf_.mesh,
                                            self._vf_.func,
                                            ftype=self._vf_.ftype,
                                            valid_time=self._vf_.valid_time,
                                            name=self._vf_.standard_properties.name
                                            )
            # this is very important as it decoupled the norm component and the vector. MUST DO THIS!

            trace_element_wise_func = dict()
            for i in safe_copy.mesh.trace.elements: # the local trace element #i on mesh boundaries
                trace_element_wise_func[i] = ___VECTOR_T_PERP_COMPONENT___(safe_copy, i)

            return _3dCSCG_VectorField(safe_copy.mesh, trace_element_wise_func,
                                       ftype='trace-element-wise',
                                       valid_time=safe_copy.valid_time,
                                       name='T-perp-component-of-' + safe_copy.standard_properties.name
                                       )

        else:
            raise NotImplementedError(f"`T_perp` component of a vector of "
                                      f"type = {self._vf_.ftype} is not implemented")



class ___VECTOR_NORM_COMPONENT___(object):
    """Here we wrap the reconstruction of standard vector such that it works like a function. Then we
    can use it to, for example, build more vectors."""
    def __init__(self, vf, i):
        """

        :param vf: the vector
        :param i: for trace element #i
        """
        self._vf_ = vf
        self._i_ = i
        self._te_ = vf.mesh.trace.elements[i]

    def __call__(self, t, xi, et, sg):
        """This actually includes a reconstruction. So we will call from xi, et, sg.

        :param t: at time t
        :param xi: must be a 1d array.
        :param et: must be a 1d array.
        :param sg: must be a 1d array.
        :return:
        """
        vf = self._vf_
        i = self._i_

        if vf.ftype == 'standard':

            xyz, w = vf.reconstruct(xi, et, sg,
                                    time=t,
                                    i=i,
                                    ravel=False,
                                    where='trace-element',
                                    structured=True)
            xyz = xyz[i]
            w = w[i]
            n = self._te_.coordinate_transformation.unit_normal_vector(xi, et, sg, parse_3_1d_eps=True)
            w_dot_n = w[0]*n[0] + w[1]*n[1] + w[2]*n[2]
            norm_component = (w_dot_n * n[0], w_dot_n * n[1], w_dot_n * n[2])
            return xyz, norm_component

        else:
            raise Exception(f"We cannot compute norm component for {self._vf_.ftype} _3dCSCG_VectorField.")

class ___VECTOR_T_PARA_COMPONENT___(object): # Trace parallel
    """Here we wrap the reconstruction of standard vector such that it works like a function. Then we
    can use it to, for example, build more vectors.

    For parallel trace component.
    """
    def __init__(self, vf, i):
        """

        :param vf: the vector
        :param i: for trace element #i
        """
        self._vf_ = vf
        self._i_ = i
        self._te_ = vf.mesh.trace.elements[i]

    def __call__(self, t, xi, et, sg):
        """This actually includes a reconstruction. So we will call from xi, et, sg.

        :param t: at time t
        :param xi: must be a 1d array.
        :param et: must be a 1d array.
        :param sg: must be a 1d array.
        :return:
        """
        vf = self._vf_
        i = self._i_

        if vf.ftype == 'standard':

            xyz, w = vf.reconstruct(xi, et, sg,
                                    time=t,
                                    i=i,
                                    ravel=False,
                                    where='trace-element',
                                    structured=True)
            xyz = xyz[i]
            w = w[i]
            n = self._te_.coordinate_transformation.unit_normal_vector(xi, et, sg, parse_3_1d_eps=True)

            w1, w2, w3 = w
            n1, n2, n3 = n

            wXn1 = w2 * n3 - w3 * n2
            wXn2 = w3 * n1 - w1 * n3
            wXn3 = w1 * n2 - w2 * n1

            n_X_lb_wXn_rb = (n2 * wXn3 - n3 * wXn2,
                             n3 * wXn1 - n1 * wXn3,
                             n1 * wXn2 - n2 * wXn1)

            return xyz, n_X_lb_wXn_rb

        else:
            raise Exception(f"We cannot compute T_para component for {self._vf_.ftype} _3dCSCG_VectorField.")

class ___VECTOR_T_PERP_COMPONENT___(object): # Trace Perpendicular
    """Here we wrap the reconstruction of standard vector such that it works like a function. Then we
    can use it to, for example, build more vectors.

    For perpendicular trace component.
    """
    def __init__(self, vf, i):
        """

        :param vf: the vector
        :param i: for trace element #i
        """
        self._vf_ = vf
        self._i_ = i
        self._te_ = vf.mesh.trace.elements[i]

    def __call__(self, t, xi, et, sg):
        """This actually includes a reconstruction. So we will call from xi, et, sg.

        :param t: at time t
        :param xi: must be a 1d array.
        :param et: must be a 1d array.
        :param sg: must be a 1d array.
        :return:
        """
        vf = self._vf_
        i = self._i_

        if vf.ftype == 'standard':

            xyz, w = vf.reconstruct(xi, et, sg,
                                    time=t,
                                    i=i,
                                    ravel=False,
                                    where='trace-element',
                                    structured=True)
            xyz = xyz[i]
            w = w[i]
            n = self._te_.coordinate_transformation.unit_normal_vector(xi, et, sg, parse_3_1d_eps=True)
            w1, w2, w3 = w
            n1, n2, n3 = n

            wXn = (w2 * n3 - w3 * n2,
                   w3 * n1 - w1 * n3,
                   w1 * n2 - w2 * n1)

            return xyz, wXn

        else:
            raise Exception(f"We cannot compute T_perp component for {self._vf_.ftype} _3dCSCG_VectorField.")












if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\field\vector.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([1,1,1], show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)], show_info=True)
    FC = FormCaller(mesh, space)

    # def p(t, x, y, z): return t + np.cos(np.pi*x) * np.cos(2*np.pi*y) * np.cos(3*np.pi*z)
    # SS = FC('scalar', p)
    # BS = FC('scalar', {'North': p, 'West':p})


    def velocity_x(t, x, y, z): return t + np.cos(2*np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z)
    def velocity_y(t, x, y, z): return t + np.sin(np.pi*x) * np.cos(2*np.pi*y) * np.sin(np.pi*z)
    def velocity_z(t, x, y, z): return t + np.sin(np.pi*x) * np.sin(np.pi*y) * np.cos(2*np.pi*z)

    SV = FC('vector', [velocity_x, velocity_y, velocity_z])
    #
    # bcDv = {'North': [velocity_x, velocity_y, velocity_z],
    #         'West': [velocity_x, velocity_y, velocity_z]}
    # BV = FC('vector', bcDv)
    # #
    # GT = SV.numerical.gradient
    # CV = SV.numerical.curl
    # DS = SV.numerical.divergence

    # print(GT)

    # flux = SV.flux
    #
    # # xi = np.linspace(-1,1,10)
    # # et = np.linspace(-1,1,10)
    # # sg = np.linspace(-1,1,10)
    #
    # flux.current_time=0
    # flux.visualize()


    norm = SV.components.norm
    norm.current_time=0
    norm.visualize()

    para = SV.components.T_para
    para.current_time=0
    para.visualize()

    perp = SV.components.T_perp
    perp.current_time=0
    perp.visualize()

    # xyz, nc = ___VECTOR_NORM_COMPONENT___(SV, 0)(0, [-1,1], [-1,0,0.5,1], [-1,0,1])
    #
    # print(np.shape(xyz), np.shape(nc))

    # print(norm)

    # print(SV.covered_mesh_boundaries)

    # GT = BV.numerical.gradient

    # BV.current_time = 1
    #
    # BV.visualize()

    # f2 = FC('2-t')
    # f2.TW.BC.body = BV
    # f2.TW.DO.push_BC_to_instant(0)

    # f2.BC.valid_boundaries=['North','West', ]
    # f2pc = f2.BC.partial_cochain
    #
    # for i in f2pc:
    #     print(i, f2pc[i])
