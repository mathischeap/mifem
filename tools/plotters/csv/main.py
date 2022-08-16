# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/16/2022 1:26 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import pandas as pd


class csvPlotter(FrozenOnly):
    """"""

    def __init__(self, csv_filename):
        """"""
        self._pdf_ = pd.read_csv(csv_filename, index_col=0)
        self._array_ = self._pdf_.to_numpy()
        self._freeze_self_()

    @property
    def pdf(self):
        return self._pdf_
    
    @property
    def array(self):
        return self._array_

    @property
    def columns(self):
        return list(self._pdf_.columns)


    
if __name__ == '__main__':
    # mpiexec -n 4 python tools/plotters/csv/main.py
    pass
