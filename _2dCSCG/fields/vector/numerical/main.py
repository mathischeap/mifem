"""

"""

from screws.freeze.main import FrozenOnly
from importlib import import_module
from screws.numerical.time_plus_2d_space.partial_derivative_as_functions import \
    NumericalPartialDerivative_txy_Functions
from _2dCSCG.fields.vector.numerical.helpers.curl import ___VECTOR_CURL_HELPER___
from _2dCSCG.fields.vector.numerical.helpers.divergence import ___VECTOR_DIV_HELPER___





class _2dCSCG_VectorField_Numerical(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    @property
    def time_derivative(self):
        """Return a _3dCSCG_ScalarField instances which is the numerical time derivative of self."""
        if self._vf_.ftype == 'standard':
            NPD4Fx = NumericalPartialDerivative_txy_Functions(self._vf_.func[0])
            NPD4Fy = NumericalPartialDerivative_txy_Functions(self._vf_.func[1])

            vector_class = getattr(import_module('_2dCSCG.fields.vector.main'), '_2dCSCG_VectorField')

            TDS = vector_class(self._vf_.mesh,
                               [NPD4Fx('t'), NPD4Fy('t')],
                               ftype='standard',
                               valid_time=self._vf_.valid_time,
                               name = 'time-derivative-of-' + self._vf_.standard_properties.name
                               )
            return TDS

        else:
            raise NotImplementedError(f"Numerical time derivative not implemented for vector type = {self._vf_.ftype}.")





    @property
    def curl(self):
        """Return a _2dCSCG_ScalarField instance which is the numerical curl of self."""
        if self._vf_.ftype == 'standard':
            u, v = self._vf_.func
            NPD3F_u = NumericalPartialDerivative_txy_Functions(u)
            NPD3F_v = NumericalPartialDerivative_txy_Functions(v)
            u_y = NPD3F_u('y')
            v_x = NPD3F_v('x')

            curl_vector = ___VECTOR_CURL_HELPER___(v_x, u_y)

            scalar_class = getattr(import_module('_2dCSCG.fields.scalar.main'), '_2dCSCG_ScalarField')

            curl_vector = scalar_class(self._vf_.mesh,
                                     curl_vector,
                                     ftype='standard',
                                     valid_time=self._vf_.valid_time,
                                     name = 'curl-of-' + self._vf_.standard_properties.name
                                     )
            return curl_vector
        else:
            raise NotImplementedError(f"Numerical curl not implemented for vector type = {self._vf_.ftype}.")




    @property
    def div(self):
        """Return a _2dCSCG_ScalarField instance which is the numerical div of self."""
        if self._vf_.ftype == 'standard':
            u, v = self._vf_.func
            NPD3F_u = NumericalPartialDerivative_txy_Functions(u)
            NPD3F_v = NumericalPartialDerivative_txy_Functions(v)
            u_x = NPD3F_u('x')
            v_y = NPD3F_v('y')

            curl_vector = ___VECTOR_DIV_HELPER___(u_x, v_y)

            scalar_class = getattr(import_module('_2dCSCG.fields.scalar.main'), '_2dCSCG_ScalarField')

            curl_vector = scalar_class(self._vf_.mesh,
                                     curl_vector,
                                     ftype='standard',
                                     valid_time=self._vf_.valid_time,
                                     name = 'div-of-' + self._vf_.standard_properties.name
                                     )
            return curl_vector
        else:
            raise NotImplementedError(f"Numerical div not implemented for vector type = {self._vf_.ftype}.")