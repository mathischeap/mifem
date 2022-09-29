# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/04 10:49 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenClass


class miUsGrid_MeshBase(FrozenClass):
    """A base for mimetic triangular or tetrahedral mesh."""

    def __init__(self, ndim, name):
        """

        Parameters
        ----------
        ndim : int
            {2 ,3}.
        """
        if ndim == 2:
            assert self.__class__.__name__ == "miUsGrid_TriangularMesh"
        elif ndim == 3:
            assert self.__class__.__name__ == "miUsGrid_TetrahedralMesh"
        else:
            raise Exception

        self.standard_properties.name = name

        self._ndim_ = ndim
        self._elements_ = None
        self._visualize_ = None
        self._boundaries_ = None
        self._freeze_self_()

    def __repr__(self):
        return f"{self.__class__.__name__}@{id(self)}"

    def __eq__(self, other):
        """"""
        return self is other

    @property
    def name(self):
        return self.standard_properties.name

    @property
    def ndim(self):
        return self._ndim_

    @property
    def elements(self):
        return self._elements_

    @property
    def visualize(self):
        return self._visualize_

    @property
    def boundaries(self):
        return self._boundaries_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
