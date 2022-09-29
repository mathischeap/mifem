# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np


class miUsGrid_SF_CochainBase(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._local_ = None
        self._freeze_self_()

    def RESET_cache(self):
        """"""

    def __getitem__(self, item):
        return self.local[item]

    def __contains__(self, item):
        return item in self.local

    def __iter__(self):
        for i in self.local:
            yield i

    def __len__(self):
        return len(self.local)


    # ------------- DEPENDENT PROPERTIES (MAJOR): When set, clear BRANCHES by set _branches_ to None -------------------
    @property
    def local(self):
        """The local cochain. Must be full. So all local mesh elements must have their local cochains!

        :return: A dict whose keys are local element indices and values are cochains (1-d arrays) in corresponding
            elements.
        :rtype: Dict[int, numpy.ndarray]
        """
        return self._local_

    @local.setter
    def local(self, local):
        numOfElements = self._sf_.mesh.elements.num.cells
        numOfBasis = self._sf_.num.basis
        try:
            assert isinstance(local, dict), \
                f"local cochain needs to be a dict, now it is a {local.__class__.__name__}."
            assert len(local) == numOfElements, \
                "local cochain has to contain cochains for all local mesh elements."
            for i, j in zip(self._sf_.mesh.elements.indices, local):
                assert np.shape(local[i]) == (numOfBasis,), \
                    f"local[{i}] shape = {np.shape(local[i])} wrong. " \
                    f"It needs to be {(numOfBasis,)}."
                assert i == j, f"mesh element index sequence is wrong."

        except AssertionError:
            raise Exception("Cannot set local cochain.")

        self.RESET_cache()
        self._local_ = local


    #--------------- DEPENDENT PROPERTIES (BRANCHES, must have the two switching methods): when set, update local ------


    #=================================== ABOVE =========================================================================

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
