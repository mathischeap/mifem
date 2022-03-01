# -*- coding: utf-8 -*-
"""3D functions."""
import numpy as np
from abc import ABC


class ScalingFunc(ABC):
    """Scaling a function: new_func = C*func."""
    def __init__(self, C):
        self._C_ = C

    def _scaled_func_(self, x, y, z):
        return self._C_ * self._func_(x, y, z)

    def __call__(self, func):
        self._func_ = func
        return self._scaled_func_






class Opposite(ABC):
    """Equal to ``ScalingFunc(-1, func)()``."""
    def __init__(self, func):
        """ """
        self._func_ = func

    def _opposite_func_(self, x, y, z):
        return -self._func_(x, y, z)

    def __call__(self):
        return self._opposite_func_






class CFG(ABC):
    """Constant function generator."""
    def __init__(self, C):
        """ """
        self._C_ = C

    def _constant_func_(self, x, y, z):
        assert np.shape(x) == np.shape(y) == np.shape(z)
        return self._C_ + np.zeros(np.shape(x))

    def __call__(self):
        return self._constant_func_







class CartSphSwitcher(ABC):
    """A spherical <-> Cartesian coordinate switcher."""
    @classmethod
    def cart2sph(cls, x, y, z):
        hxy = np.hypot(x, y)
        r = np.hypot(hxy, z)
        el = np.arctan2(z, hxy)
        az = np.arctan2(y, x)
        return az, el, r

    @classmethod
    def sph2cart(cls, az, el, r):
        r_cos_theta = r * np.cos(el)
        x = r_cos_theta * np.cos(az)
        y = r_cos_theta * np.sin(az)
        z = r * np.sin(el)
        return x, y, z


class CartCylSwitcher(ABC):
    """A cylinder <-> Cartesian coordinate switcher."""
    @classmethod
    def cart2cyl(cls, x, y, z):
        rho = np.sqrt(x ** 2 + y ** 2)
        phi = np.arctan2(y, x)
        return rho, phi, z

    @classmethod
    def cyl2cart(cls, rho, phi, z):
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return x, y, z



def _0_(x, y, z):
    assert np.shape(x) == np.shape(y) == np.shape(z)
    return np.zeros(np.shape(x))





def angle_between_two_vectors(v1, v2):
    """Compute the angle between the angle between v1 and v2 (can be
    of any dimension).

    :param v1:
    :param v2:
    :type v1: tuple, list, np.array
    :type v2: tuple, list, np.array
    :return: a float between [0, pi]
    """
    dot_product = sum([a*b for a, b in zip(v1, v2)])
    len_v1 = np.sqrt(np.sum([a**2 for a in v1]))
    len_v2 = np.sqrt(np.sum([a**2 for a in v2]))

    _ = dot_product / (len_v1 * len_v2)

    if _ > 1:
        _ = 1
    elif _ < -1:
        _ = -1
    else:
        pass

    return np.arccos(_)