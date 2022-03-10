import numpy as np
from screws.freeze.main import FrozenOnly
from screws.decorators.accepts import accepts
from screws.numerical._2d_space.Jacobian_21 import NumericalJacobian_xy_22


class InterpolationBase(FrozenOnly):
    """ """
    @accepts('self', 'Region')
    def __init__(self, region):
        self._region_ = region
        self._ndim_ = region.ndim

    @property
    def ndim(self):
        return self._ndim_

    @staticmethod
    def ___check_rs___(r, s):
        """ r, s be in [0, 1]. """
        if r.__class__.__name__ == 'ndarray':
            pass
        elif isinstance(r, (int, float)):
            r = np.array([r])
        elif isinstance(r, list):
            r = np.array(r)
        else:
            raise Exception()

        if s.__class__.__name__ == 'ndarray':
            pass
        elif isinstance(s, (int, float)):
            s = np.array([s])
        elif isinstance(s, list):
            s = np.array(s)
        else:
            raise Exception()

        assert np.shape(r) == np.shape(s), " <Interpolation> : inputs shape dismatch."
        return r, s

    def mapping(self, r, s):
        raise NotImplementedError()

    def mapping_X(self, r, s):
        raise NotImplementedError()
    def mapping_Y(self, r, s):
        raise NotImplementedError()


    def __call__(self, r, s):
        return self.mapping(r, s)


    def Jacobian_matrix(self, r, s):
        """
        r, s, t be in [0, 1].

        This is a general Jacobian_matrix using numerical derivative. To avoid
        this, just overwrite it in the child class.

        """
        r, s = self.___check_rs___(r, s)
        NJ22 = NumericalJacobian_xy_22(self.mapping)
        return NJ22.scipy_derivative(r, s)

    def Jacobian_Xr(self, r, s):
        raise NotImplementedError()
    def Jacobian_Xs(self, r, s):
        raise NotImplementedError()

    def Jacobian_Yr(self, r, s):
        raise NotImplementedError()
    def Jacobian_Ys(self, r, s):
        raise NotImplementedError()

