# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 1:06 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly




class mpRfT2_CellTypeWr2Metric_Base(FrozenOnly):
    """
    A base for all cell types w.r.t. metric. For each type of cscg mesh element, we can classify its
    sub-cells into different types. These types are all coded here.
    """
    @property
    def mark(self):
        """
        A mark is key that identifies the cell metric. If the marks of two cells are the same, then
        they have the same metric, otherwise, their metric are different. A mark normally is a string. But for
        chaotic cell, it is an int: the id of the object.

        :return:
        """
        # noinspection PyUnresolvedReferences
        return self._mark_

    def __eq__(self, other):
        """We ask that the marks to be the same."""
        # The later judge is to make sure we are not comparing it to something else having a mark property
        return self.mark == other.mark and other.___IS_mpRfT2_CellTypeWr2Metric___

    @property
    def ___IS_mpRfT2_CellTypeWr2Metric___(self):
        return True




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
