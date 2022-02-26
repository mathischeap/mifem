# -*- coding: utf-8 -*-
"""2D numerical."""
from types import FunctionType, MethodType
import numpy as np
from abc import ABC
from scipy.misc import derivative
from screws.numerical._1d import NumericalDerivative_fx






class NumericalJacobian_xy_t_21(ABC):
    """For a mapping: ``XY(t) = (x, y) = (X(t), Y(t))``, We compute ``dx/dt``, and ``dy/dt``.
    """
    def __init__(self, func21):
        """ """
        self._func21_ = func21

    def ___evaluate_func21_for_x_t___(self, t):
        return self._func21_(t)[0]
    def ___evaluate_func21_for_y_t___(self, t):
        return self._func21_(t)[1]

    def scipy_derivative(self, t, dt=1e-6, n=1, order=3):
        Xt = NumericalDerivative_fx(self.___evaluate_func21_for_x_t___, t,
                                    dx=dt, n=n, order=order).scipy_derivative()
        Yt = NumericalDerivative_fx(self.___evaluate_func21_for_y_t___, t,
                                    dx=dt, n=n, order=order).scipy_derivative()
        return Xt, Yt

    def check_Jacobian(self, Jacobian, t, tolerance=1e-6):
        """Check if ``Jacobian(t) == self.scipy_derivative(t)`` at nodes ``t``. """
        self_J = self.scipy_derivative(t)
        give_J = Jacobian(t)
        result = [None, None]
        for i in range(2):
            absolute_error = np.max(np.abs(self_J[i]-give_J[i]))
            if absolute_error < tolerance:
                result[i] = True
            else:
                relative_error = np.max(np.abs((self_J[i]-give_J[i])/self_J[i]))
                if relative_error < tolerance:
                    result[i] = True
                else:
                    result[i] = False
        return tuple(result)






class NumericalJacobian_xy_22(ABC):
    """
    For a mapping: ``x = Phi_x(r, s), y = Phi_y(r, s)``,
    ``self._func_(r, s) = (Phi_x(r, s), Phi_y(r, s))``, we compute the its Jacobian numerically:
    ``(( dx/dr, dx/ds ), ( dy/dr, dy/ds ))``.

    """
    def __init__(self, func22):
        """ """
        self._func22_ = func22

    def ___evaluate_func22_for_x_rs___(self, r, s):
        return self._func22_(r, s)[0]
    def ___evaluate_func22_for_y_rs___(self, r, s):
        return self._func22_(r, s)[1]

    def scipy_derivative(self, r, s, dr_ds=1e-8, n=1, order=3):
        xr, xs = NumericalPartialDerivative_xy(self.___evaluate_func22_for_x_rs___,
                                               r, s, dx_dy=dr_ds, n=n, order=order).scipy_total
        yr, ys = NumericalPartialDerivative_xy(self.___evaluate_func22_for_y_rs___,
                                               r, s, dx_dy=dr_ds, n=n, order=order).scipy_total
        return ((xr, xs),
                (yr, ys))





class NumericalPartialDerivative_xy(ABC):
    """
    Numerical partial derivative; we call it '2' because we compute a function or method that
    like: ``a=f(x,y)``.

    :param func:
    :param x:
    :param y:
    :param dx_dy: The interval. The smaller, the more accurate.
    :param n: nth order derivative.
    :param order: How many points are used to approximate the derivative.
    """
    def __init__(self, func, x, y, dx_dy=1e-6, n=1, order=3):
        self.___check_func___(func)
        self.___check_xy___(x, y)
        self.___check_dx_dy___(dx_dy)
        self.___check_n___(n)
        self.___check_order___(order)

    def ___check_func___(self, func):
        """ """
        assert callable(func), " <PartialDerivative> : func is not callable."
        if isinstance(func, FunctionType):
            # noinspection PyUnresolvedReferences
            assert func.__code__.co_argcount == 2, " <PartialDerivative> : need a func of 2 args."
        elif isinstance(func, MethodType):
            # noinspection PyUnresolvedReferences
            assert func.__code__.co_argcount == 3, \
                " <PartialDerivative> : need a method of 2 args (3 including self)."
        elif func.__class__.__name__ == 'partial':
            # noinspection PyUnresolvedReferences
            if isinstance(func.func, FunctionType):
                # noinspection PyUnresolvedReferences
                assert func.func.__code__.co_argcount == 3
            elif isinstance(func.func, MethodType):
                # noinspection PyUnresolvedReferences
                assert func.func.__code__.co_argcount == 4
            else:
                raise Exception()
        else:
            raise NotImplementedError(func.__class__.__name__)
        self._func_ = func

    def ___check_xy___(self, x, y):
        """ """
        self._x_, self._y_ = x, y
        assert np.shape(self._x_) == np.shape(self._y_), \
            " <PartialDerivative> : xy of different shapes."

    def ___check_dx_dy___(self, dx_dy):
        """ """
        if isinstance(dx_dy, (int, float)):
            self._dx_ = self._dy_ = dx_dy
        else:
            assert np.shape(dx_dy) == (2,), " <PartialDerivative> : dx_dy shape wrong."
            self._dx_, self._dy_ = dx_dy
        assert all([isinstance(d, (int, float)) and d > 0 for d in (self._dx_, self._dy_)])

    def ___check_n___(self, n):
        """ """
        assert n % 1 == 0 and n >= 1, " <PartialDerivative> : n = {} is wrong.".format(n)
        self._n_ = n

    def ___check_order___(self, order):
        """ """
        assert order % 2 == 1 and order > 0, " <PartialDerivative> : order needs to be odd positive."
        self._order_ = order

    def ___evaluate_func_for_x___(self, x):
        return self._func_(x, self._y_)

    def ___evaluate_func_for_y___(self, y):
        return self._func_(self._x_, y)

    def scipy_partial(self, d_):
        """ We compute ``df/d_`` at points ``*xyz.``"""
        if d_ == 'x':
            # noinspection PyTypeChecker
            return derivative(self.___evaluate_func_for_x___, self._x_, dx=self._dx_,
                              n=self._n_, order=self._order_)
        elif d_ == 'y':
            # noinspection PyTypeChecker
            return derivative(self.___evaluate_func_for_y___, self._y_, dx=self._dy_,
                              n=self._n_, order=self._order_)
        else:
            raise Exception(" <PartialDerivative> : dx or dy or dz? ")

    @property
    def scipy_total(self):
        px = self.scipy_partial('x')
        py = self.scipy_partial('y')
        return px, py

    def check_partial_x(self, px_func, tolerance=1e-5):
        self_px = self.scipy_partial('x')
        func_px = px_func(self._x_, self._y_)
        absolute_error = np.max(np.abs(func_px-self_px))
        if absolute_error < tolerance:
            return True
        relative_error = np.max(np.abs((func_px-self_px)/self_px))
        if relative_error < tolerance:
            return True
        else:
            return False

    def check_partial_y(self, py_func, tolerance=1e-5):
        self_py = self.scipy_partial('y')
        func_py = py_func(self._x_, self._y_)
        absolute_error = np.max(np.abs(func_py-self_py))
        if absolute_error < tolerance:
            return True
        relative_error = np.max(np.abs((func_py-self_py)/self_py))
        if relative_error < tolerance:
            return True
        else:
            return False

    def check_total(self, px_func, py_func, tolerance=1e-5):
        return (self.check_partial_x(px_func, tolerance=tolerance),
                self.check_partial_y(py_func, tolerance=tolerance))