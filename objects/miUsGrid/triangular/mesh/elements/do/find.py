# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/19 9:52 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import SIZE, np


class miUsGrid_TriangularMesh_Elements_DO_FIND(FrozenOnly):
    """"""

    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()


    def rank_of_element(self, i: int) -> int:
        DISTRI = self._elements_.distributions
        if isinstance(i, str): i = int(i)
        if SIZE <= 4:
            for nC in range(SIZE):
                if i in DISTRI[nC]: return nC
            raise Exception()

        midCore0 = 0
        midCore1 = SIZE // 2
        midCore2 = SIZE
        while i not in DISTRI[midCore1] and midCore1 - midCore0 > 2 and midCore2 - midCore1 > 2:
            if i > max(DISTRI[midCore1]):
                midCore0 = midCore1
                midCore1 = (midCore0 + midCore2) // 2
            elif i < min(DISTRI[midCore1]):
                midCore2 = midCore1
                midCore1 = (midCore0 + midCore2) // 2
            else:
                raise Exception
        if i in DISTRI[midCore1]:
            return midCore1
        elif i > np.max(DISTRI[midCore1]):
            for noCore in range(midCore1, midCore2):
                if i in DISTRI[noCore]: return noCore
        elif i < np.min(DISTRI[midCore1]):
            for noCore in range(midCore0, midCore1):
                if i in DISTRI[noCore]: return noCore
        else:
            raise Exception



if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/elements/do/find.py
    pass