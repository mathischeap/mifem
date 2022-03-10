
from screws.freeze.main import FrozenOnly
from importlib import import_module





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

            vector_class = getattr(import_module('_3dCSCG.fields.vector.main'), '_3dCSCG_VectorField')

            safe_copy = vector_class(self._vf_.mesh,
                                            self._vf_.func,
                                            ftype=self._vf_.ftype,
                                            valid_time=self._vf_.valid_time,
                                            name=self._vf_.standard_properties.name
                                            )
            # this is very important as it decoupled the norm component and the vector. MUST do THIS!

            trace_element_wise_func = dict()
            for i in safe_copy.mesh.trace.elements: # the local trace element #i on mesh boundaries
                trace_element_wise_func[i] = ___VECTOR_NORM_COMPONENT___(safe_copy, i)

            return vector_class(safe_copy.mesh, trace_element_wise_func,
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

            vector_class = getattr(import_module('_3dCSCG.fields.vector.main'), '_3dCSCG_VectorField')

            safe_copy = vector_class(self._vf_.mesh,
                                            self._vf_.func,
                                            ftype=self._vf_.ftype,
                                            valid_time=self._vf_.valid_time,
                                            name=self._vf_.standard_properties.name
                                            )
            # this is very important as it decoupled the norm component and the vector. MUST do THIS!

            trace_element_wise_func = dict()
            for i in safe_copy.mesh.trace.elements: # the local trace element #i on mesh boundaries
                trace_element_wise_func[i] = ___VECTOR_T_PARA_COMPONENT___(safe_copy, i)

            return vector_class(safe_copy.mesh, trace_element_wise_func,
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

            vector_class = getattr(import_module('_3dCSCG.fields.vector.main'), '_3dCSCG_VectorField')

            safe_copy = vector_class(self._vf_.mesh,
                                            self._vf_.func,
                                            ftype=self._vf_.ftype,
                                            valid_time=self._vf_.valid_time,
                                            name=self._vf_.standard_properties.name
                                            )
            # this is very important as it decoupled the norm component and the vector. MUST do THIS!

            trace_element_wise_func = dict()
            for i in safe_copy.mesh.trace.elements: # the local trace element #i on mesh boundaries
                trace_element_wise_func[i] = ___VECTOR_T_PERP_COMPONENT___(safe_copy, i)

            return vector_class(safe_copy.mesh, trace_element_wise_func,
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

            # #------- option 1 ----------------------------------------------------------------------
            wXn = (w2 * n3 - w3 * n2,
                   w3 * n1 - w1 * n3,
                   w1 * n2 - w2 * n1)
            #--------- option 2: TEST FOR ORTHOGONAL MESH ONLY ------------------------------------
            # side = self._te_.CHARACTERISTIC_side
            # if side in 'NS':
            #     wXn =  (w1, w3, -w2)
            # elif side in 'WE':
            #     wXn =  (-w3, w2, w1)
            # elif side in 'BF':
            #     wXn =  (w2, -w1, w3)
            # else:
            #     raise Exception

            #=========================================================================================


            return xyz, wXn

        else:
            raise Exception(f"We cannot compute T_perp component for {self._vf_.ftype} _3dCSCG_VectorField.")