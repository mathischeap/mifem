# -*- coding: utf-8 -*-
"""2D functions."""
import numpy as np
from abc import ABC


class ScalingFunc(ABC):
    """
    Scaling a function: new_func = C*func.

    To scale a function, saying `func`, with constant C, do for example:

    .. doctest::

        >>> new_func = ScalingFunc(3)(CFG(5)())
        >>> new_func(0, 0)
        15.0

    ``newfunc`` is the scaled function which alway return 15.
    """

    def __init__(self, C):
        self._C_ = C

    def _scaled_func_(self, x, y):
        return self._C_ * self._func_(x, y)

    def __call__(self, func):
        self._func_ = func
        return self._scaled_func_


class Opposite(ABC):
    """
    Equal to ``ScalingFunc(-1)(func)``.

    .. doctest::

        >>> new_func = Opposite(CFG(5)())()
        >>> new_func(0, 0)
        -5.0
    """

    def __init__(self, func):
        self._func_ = func

    def _opposite_func_(self, x, y):
        return -self._func_(x, y)

    def __call__(self):
        return self._opposite_func_


def _0_(x, y):
    """
    A function always returns ``0``.

    :param x:
    :param y:
    :return:

    .. doctest::

        >>> _0_(100,1000)
        array(0.)
    """
    assert np.shape(x) == np.shape(y)
    return np.zeros(np.shape(x))

def _0t_(t, x, y):
    """
    A function always returns ``0``.

    :param x:
    :param y:
    :return:

    .. doctest::

        >>> _0t_(100, 100, 1000)
        0.0
    """
    assert np.shape(x) == np.shape(y)
    return np.zeros(np.shape(x)) + 0 * t


class CFG(ABC):
    """
    Constant function generator: ``CF = CFG(5)()``.

    .. doctest::

        >>> cf = CFG(10.5)()
        >>> cf(1,1)
        10.5
    """

    def __init__(self, C):
        """ """
        self._C_ = C

    def _constant_func_(self, x, y):
        assert np.shape(x) == np.shape(y)
        return self._C_ + np.zeros(np.shape(x))

    def __call__(self):
        return self._constant_func_


class CFGt(ABC):
    """
    Constant function generator: ``CF = CFG(5)()``.

    .. doctest::

        >>> cf = CFGt(10.5)()
        >>> cf(1, 1,1)
        10.5
    """

    def __init__(self, C):
        """ """
        self._C_ = C

    def _constant_func_(self, t, x, y):
        assert np.shape(x) == np.shape(y)
        return self._C_ + np.zeros(np.shape(x)) + 0 * t

    def __call__(self):
        return self._constant_func_






class CartPolSwitcher(ABC):
    """
    A polar <-> Cartesian coordinates switcher.

    .. doctest::

        >>> CartPolSwitcher.cart2pol(2,2)
        (2.8284271247461903, 0.7853981633974483)
        >>> CartPolSwitcher.pol2cart(1, np.pi/4)
        (0.7071067811865476, 0.7071067811865476)
    """

    @classmethod
    def cart2pol(cls, x, y):
        rho = np.sqrt(x ** 2 + y ** 2)
        phi = np.arctan2(y, x)
        return rho, phi

    @classmethod
    def pol2cart(cls, rho, phi):
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return x, y




def angle(origin, pt):
    """
    Angle between the vector from origin to pt and the x-direction vector.

    For example:

    .. doctest::

        >>> angle((0,0), (1,1)) # should return pi/4
        0.7853981633974484
    """
    x1, y1 = (1, 0)
    x2, y2 = (pt[0] - origin[0], pt[1] - origin[1])
    inner_product = x1 * x2 + y1 * y2
    len1 = np.hypot(x1, y1)
    len2 = np.hypot(x2, y2)
    if y2 < 0:
        return 2 * np.pi - np.arccos(inner_product / (len1 * len2))
    else:
        return np.arccos(inner_product / (len1 * len2))


def distance(p1, p2):
    """ Compute distance between two points. """
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def if_two_lines_parallel(a1, a2, b1, b2):
    """
    Check if line p1-p2 is parallel with line p3-p4. If they are parallel but pointing different direction,
    return False

    :param a1:
    :param a2:
    :param b1:
    :param b2:
    :return:
    :rtype: bool
    """
    angle1 = angle(a1, a2)
    angle2 = angle(b1, b2)
    return angle1 == angle2


def __det__(a, b):
    return a[0] * b[1] - a[1] * b[0]

def find_line_intersection(a1, a2, b1, b2):
    """

    :param a1:
    :param a2:
    :param b1:
    :param b2:
    :return:

    .. doctest::

        >>> a1, a2, b1, b2 = (0,0), (1,0), (0,1), (1,1)
        >>> find_line_intersection(a1, a2, b1, b2)
        >>> a1, a2, b1, b2 = (0,0), (1,1), (1,0), (0,1)
        >>> find_line_intersection(a1, a2, b1, b2)
        (0.5, 0.5)
        >>> a1, a2, b1, b2 = (0,0), (1,0), (0,1), (1,2)
        >>> find_line_intersection(a1, a2, b1, b2)
        (-1.0, 0.0)
        >>> find_line_intersection(a1, a2, b2, b1)
        (-1.0, -0.0)
    """
    line1, line2 = (a1, a2), (b1, b2)
    xdiff = (a1[0] - a2[0], b1[0] - b2[0])
    ydiff = (a1[1] - a2[1], b1[1] - b2[1])

    div = __det__(xdiff, ydiff)
    if div == 0:
       return None

    d = (__det__(*line1), __det__(*line2))
    x = __det__(d, xdiff) / div
    y = __det__(d, ydiff) / div
    return x, y







if __name__ == "__main__":
    import doctest
    doctest.testmod()