# -*- coding: utf-8 -*-
import numpy as np
from screws.freeze.main import FrozenOnly
from screws.decorators.accepts import accepts
from screws.numerical._2d_space.Jacobian_21 import NumericalJacobian_xy_22
from pynverse import inversefunc

class InterpolationBase(FrozenOnly):
    """ """
    @accepts('self', 'Region')
    def __init__(self, region):
        self._region_ = region
        self._ndim_ = region.ndim
        self._Rx00_ = None
        self._Sy00_ = None

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

        assert np.shape(r) == np.shape(s), " <Interpolation> : inputs shape dis-match."
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




    def mapping_Xr_at_s0(self, r):
        """x = mapping_X(r, 0)"""
        return self.mapping_X(r, np.zeros_like(r))
    def mapping_Ys_at_r0(self, s):
        """y = mapping_Y(0, s)"""
        return self.mapping_Y(np.zeros_like(s), s)





    def ___mapping_Xr_s0___(self, r):
        """x = mapping_X(r, 0); do not use, for inverse only"""
        return self.mapping_X(r, 0)
    def ___mapping_Ys_r0___(self, s):
        """y = mapping_Y(0, s): do n ot use, for inverse only"""
        return self.mapping_Y(0, s)


    def ___inverse_mapping_r_x_s0___(self, x):
        """Return r according to the inverse function of x = mapping_X(r, 0 , 0)

        Parameters
        ----------
        x :

        Returns
        -------
        r :

        """
        if self._Rx00_ is None:
            self._Rx00_ = inversefunc(self.___mapping_Xr_s0___)
        return self._Rx00_(x)

    def ___inverse_mapping_s_y_r0___(self, y):
        """Return s according to the inverse function of y = mapping_Y(0, s , 0)

        Parameters
        ----------
        y :

        Returns
        -------
        s :

        """
        if self._Sy00_ is None:
            self._Sy00_ = inversefunc(self.___mapping_Ys_r0___)
        return self._Sy00_(y)
