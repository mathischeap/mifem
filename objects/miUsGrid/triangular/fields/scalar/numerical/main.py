# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly
from importlib import import_module
from screws.numerical.time_plus_2d_space.partial_derivative_as_functions import \
    NumericalPartialDerivative_txy_Functions

from objects.CSCG._2d.fields.scalar.numerical.helpers.curl import ___VECTOR_CURL_HELPER___



class miUsGrid_Triangular_Scalar_Numerical(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    @property
    def time_derivative(self):
        """Return a miUsGrid_Triangular_Scalar instances which is the numerical time derivative of self."""
        if self._sf_.ftype == 'standard':
            NPD4F = NumericalPartialDerivative_txy_Functions(self._sf_.func[0])

            TDS = self._sf_.__class__(self._sf_.mesh, NPD4F('t'),
                              ftype='standard',
                              valid_time=self._sf_.valid_time,
                              name = 'time-derivative-of-' + self._sf_.standard_properties.name
                              )
            return TDS

        else:
            raise NotImplementedError(f"Numerical time derivative not implemented for scalar type = {self._sf_.ftype}.")


    @property
    def grad(self):
        """Return a miUsGrid_Triangular_Vector instance which is the numerical gradient of self."""
        if self._sf_.ftype == 'standard':
            NPD3F_f = NumericalPartialDerivative_txy_Functions(self._sf_.func[0])
            f_x = NPD3F_f('x')
            f_y = NPD3F_f('y')

            base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-4]) + '.'

            vector_class = getattr(import_module(base_path + 'vector.main'), 'miUsGrid_Triangular_Vector')

            grad_vector = vector_class(self._sf_.mesh,
                                     (f_x, f_y),
                                     ftype='standard',
                                     valid_time=self._sf_.valid_time,
                                     name = 'gradient-of-' + self._sf_.standard_properties.name
                                     )
            return grad_vector

        else:
            raise NotImplementedError(f"Numerical gradient not implemented for scalar type = {self._sf_.ftype}.")


    @property
    def curl(self):
        """Return a miUsGrid_Triangular_Vector instance which is the numerical curl of self."""
        if self._sf_.ftype == 'standard':
            NPD3F_f = NumericalPartialDerivative_txy_Functions(self._sf_.func[0])
            f_x = NPD3F_f('x')
            f_y = NPD3F_f('y')
            m_f_x = ___VECTOR_CURL_HELPER___(f_x)

            base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-4]) + '.'

            vector_class = getattr(import_module(base_path + 'vector.main'), 'miUsGrid_Triangular_Vector')

            grad_vector = vector_class(self._sf_.mesh,
                                     (f_y, m_f_x),
                                     ftype='standard',
                                     valid_time=self._sf_.valid_time,
                                     name = 'curl-of-' + self._sf_.standard_properties.name
                                     )
            return grad_vector

        else:
            raise NotImplementedError(f"Numerical curl not implemented for scalar type = {self._sf_.ftype}.")