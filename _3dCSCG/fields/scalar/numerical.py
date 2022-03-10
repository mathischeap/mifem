


from screws.freeze.main import FrozenOnly
from importlib import import_module
from screws.numerical.time_plus_3d_space.partial_derivative_as_functions import NumericalPartialDerivative_txyz_Functions


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
            scalar_class = getattr(import_module('_3dCSCG.fields.scalar.main'), '_3dCSCG_ScalarField')

            TDS = scalar_class(self._sf_.mesh, NPD4F('t'),
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
            vector_class = getattr(import_module('_3dCSCG.fields.vector.main'), '_3dCSCG_VectorField')

            GV = vector_class(self._sf_.mesh, (NPD4F('x'), NPD4F('y'), NPD4F('z')),
                              ftype='standard',
                              valid_time=self._sf_.valid_time,
                              name = 'gradient-of-' + self._sf_.standard_properties.name
                              )
            return GV
        else:
            raise NotImplementedError(f"Numerical gradient not implemented for scalar type = {self._sf_.ftype}.")


