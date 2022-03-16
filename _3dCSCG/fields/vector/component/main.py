
from screws.freeze.main import FrozenOnly
from importlib import import_module

from _3dCSCG.fields.vector.component.helpers.norm import ___VECTOR_NORM_COMPONENT___
from _3dCSCG.fields.vector.component.helpers.T_para import ___VECTOR_T_PARA_COMPONENT___
from _3dCSCG.fields.vector.component.helpers.T_perp import ___VECTOR_T_PERP_COMPONENT___



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

        This is a norm component, thus it is a vector. When it is a flux, then it becomes a scalar.

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