# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/16/2022 1:26 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
import pandas as pd

from tools.filer.csv.visualize.main import csvFilerVisualize
from tools.filer.csv.do.main import csvFilerDo
from root.config.main import SIZE


class csvFiler(FrozenOnly):
    """We'd better do not use this class in parallel, although it works."""
    def __init__(self, csv_filename):
        assert SIZE == 1, f"csvFiler better not work in parallel."
        assert isinstance(csv_filename, str), "csv filename must be a str."
        if csv_filename[-4:] != '.csv': csv_filename += '.csv'
        self._df_ = pd.read_csv(csv_filename, index_col=0)
        self._visualize_ = csvFilerVisualize(self)
        self._do_ = csvFilerDo(self)
        self._freeze_self_()

    @property
    def df(self):
        return self._df_

    @property
    def columns(self):
        return self._df_.columns

    @property
    def visualize(self):
        return self._visualize_

    @property
    def do(self):
        return self._do_

if __name__ == '__main__':
    # mpiexec -n 1 python tools/filer/csv/main.py
    import os
    current_dir = os.path.dirname(__file__)
    csv = csvFiler(current_dir + '\csv_test')
    csv.do.drop(0)
    csv.df['enstrophy'] = csv.df['enstrophy'] - csv.df.loc[1,'enstrophy']

    csv.visualize.plot('t', 'enstrophy', style='-',
                        yticks = [-1e-13, 0, 1e-13,],
                        xlabel=r"$t$",
                        ylabel = r'$\mathcal{E}^h$',
                        )