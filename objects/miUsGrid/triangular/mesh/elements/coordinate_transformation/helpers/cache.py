# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/6/2022 6:05 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
import numpy as np


class miUsTriangle_Elements_CT_Cache(FrozenOnly):
    """"""

    def __init__(self, elements, method, xi ,et, intermediate_data=None):
        """

        Parameters
        ----------
        elements
        method
        xi
        et
        intermediate_data
        """
        if isinstance(xi, (int, float)): xi = np.array([xi,])
        if isinstance(et, (int, float)): et = np.array([et,])

        if xi.__class__.__name__ != 'ndarray': xi = np.array(xi)
        if et.__class__.__name__ != 'ndarray': et = np.array(et)
        assert xi.shape == et.shape, f"xi, eta shape dis-match!"

        self._elements_ = elements
        self._method_ = method
        self._xi_ = xi
        self._et_ = et
        self._cache_key_ = elements.___Pr_EWC_cache_key___
        if isinstance(self._cache_key_, str) or method == 'mapping': # no cache
            self._no_cache_ = True
        else:
            self._no_cache_ = False

        self._cache_ = dict()
        if intermediate_data is not None:
            assert intermediate_data.__class__.__name__ == "miUsTriangle_Elements_CT_Cache"
            assert intermediate_data._elements_._mesh_ == self._elements_._mesh_
            assert np.all(intermediate_data._xi_ == xi) and np.all(intermediate_data._et_ == et)
        else:
            pass

        self._itmD_ = intermediate_data
        self._freeze_self_()



    def __getitem__(self, item):
        """

        Parameters
        ----------
        item

        Returns
        -------

        """
        element = self._elements_[item]
        if self._no_cache_:
            # wo do not use intermediate data as it is not cached anyway.
            return getattr(element.coordinate_transformation, self._method_)(self._xi_, self._et_)

        else:
            cKey = self._cache_key_(item)
            if cKey in self._cache_:
                return self._cache_[cKey]
            else:
                if self._itmD_ is None:
                    data = getattr(element.coordinate_transformation, self._method_)(self._xi_, self._et_)
                else:
                    data = getattr(element.coordinate_transformation, self._method_)(
                        self._xi_, self._et_, item = self._itmD_[item]
                    )
                self._cache_[cKey] = data
                return data


if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/elements/coordinate_transformation/helpers/cache.py
    pass
