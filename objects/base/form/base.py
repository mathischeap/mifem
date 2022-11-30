# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 10:36 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenClass


class FormBase(FrozenClass):
    """A base of all forms."""

    def __init__(self, mesh, space, name):
        """"""
        self._mesh_ = mesh
        self._space_ = space
        self.standard_properties.___PRIVATE_add_tag___('form')
        self.standard_properties.name = name

        self._CF_ = None
        self._BC_ = None

    @property
    def mesh(self):
        """Return the mesh."""
        return self._mesh_

    @property
    def space(self):
        """Return the basis function space."""
        return self._space_

    @property
    def p(self):
        """Return the degree of basis functions."""
        return self.space.p

    @property
    def name(self):
        return self.standard_properties.name

    @property
    def CF(self):
        """Continuous Form."""
        return self._CF_

    @CF.setter
    def CF(self, CF):
        self.___Pr_check_CF___(CF)
        self._CF_ = CF

    def ___Pr_check_CF___(self, CF):
        raise NotImplementedError()


    @property
    def BC(self):
        return self._BC_

    def ___Pr_check_BC_CF___(self, CF):
        raise NotImplementedError()





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
