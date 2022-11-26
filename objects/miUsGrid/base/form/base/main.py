# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('../')

from components.freeze.main import FrozenClass

class miUsGrid_FormBase(FrozenClass):
    """"""

    def __init__(self, mesh, space, name):
        """"""
        self._mesh_ = mesh
        self._space_ = space
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('form')
        self.standard_properties.___PRIVATE_add_tag___('miUsGrid_form')

        self._CF_ = None
        self._discretize_ = None
        self._reconstruct_ = None
        self._cochain_ = None
        self._num_ = None

    @property
    def name(self):
        return self.standard_properties.name

    @property
    def mesh(self):
        return self._mesh_

    @property
    def space(self):
        return self._space_

    @property
    def ndim(self):
        return self.mesh.ndim

    @property
    def CF(self):
        """Continuous Form."""
        return self._CF_

    @CF.setter
    def CF(self, CF):
        assert self._CF_ is None, f"Already set CF, it cannot be changed."
        self.___Pr_check_CF___(CF)
        self._CF_ = CF

    def ___Pr_check_CF___(self, CF):
        raise NotImplementedError()


    @property
    def reconstruct(self):
        return self._reconstruct_

    @property
    def discretize(self):
        return self._discretize_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def num(self):
        return self._num_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
