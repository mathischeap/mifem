# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/16/2022 10:27 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class csvFilerDo(FrozenOnly):
    """"""

    def __init__(self, filer):
        """"""
        self._filer_ = filer
        self._freeze_self_()


    def drop(self, index):
        """delete a row indexed `index`.

        Parameters
        ----------
        index

        Returns
        -------

        """
        self._filer_._df_ = self._filer_.df.drop(index)





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
