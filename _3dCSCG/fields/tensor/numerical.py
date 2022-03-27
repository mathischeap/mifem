

from screws.freeze.main import FrozenOnly

from importlib import import_module
from screws.numerical.time_plus_3d_space.partial_derivative_as_functions import \
    NumericalPartialDerivative_txyz_Functions


class _3dCSCG_TensorField_Numerical(FrozenOnly):
    """"""
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
            tensor_class = getattr(import_module('_3dCSCG.fields.tensor.main'), '_3dCSCG_TensorField')
            TDT = tensor_class(self._tf_.mesh,
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
            vector_class = getattr(import_module('_3dCSCG.fields.vector.main'), '_3dCSCG_VectorField')
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
    """"""
    def __init__(self, fx, fy, fz):
        self._fx_ = fx
        self._fy_ = fy
        self._fz_ = fz

    def __call__(self, t, x, y, z):
        return self._fx_(t, x, y, z) + self._fy_(t, x, y, z) + self._fz_(t, x, y, z)


