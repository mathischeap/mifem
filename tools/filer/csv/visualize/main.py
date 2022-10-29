# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/16 6:03 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from screws.matplotWrappers.plot import plot, semilogy, loglog
from root.config.main import RANK, MASTER_RANK

class csvFilerVisualize(FrozenOnly):
    """"""

    def __init__(self, filer):
        """"""
        self._filer_ = filer
        self._freeze_self_()

    def plot(self, column_x,  column_y, **kwargs):
        """We plot data of 'column_x' (x-axis) and data of 'column_y' (y-axis).

        Parameters
        ----------
        column_x
        column_y : str, list, tuple
            A column or a list of columns.
        kwargs

        Returns
        -------

        """
        if RANK != MASTER_RANK: return

        if isinstance(column_x, str):
            assert column_x in self._filer_.columns, f"column_x = {column_x} is not valid."
            x = self._filer_.df[column_x]
            if isinstance(column_y, str):
                # plot one single line.
                assert column_y in self._filer_.columns, f"column_y = {column_y} is not valid."
                y = self._filer_.df[column_y]
                return plot(x, y, **kwargs)
            elif isinstance(column_y, (list, tuple)):
                X = list()
                Y = list()
                for i, col_y in enumerate(column_y):
                    assert col_y in self._filer_.columns, f"column_y[{i}] = {col_y} is not valid."
                    X.append(x)
                    Y.append(self._filer_.df[col_y])
                return plot(X, Y, num_lines=i+1, **kwargs)
            else:
                raise Exception()
        else:
            raise NotImplementedError()

    def semilogy(self, column_x,  column_y, **kwargs):
        """We plot data of 'column_x' (x-axis) and data of 'column_y' (y-axis).

        Parameters
        ----------
        column_x
        column_y
        kwargs

        Returns
        -------

        """
        if RANK != MASTER_RANK: return

        if isinstance(column_x, str):
            assert column_x in self._filer_.columns, f"column_x = {column_x} is not valid."
            x = self._filer_.df[column_x]
            if isinstance(column_y, str):
                # plot one single line.
                assert column_y in self._filer_.columns, f"column_y = {column_y} is not valid."
                y = self._filer_.df[column_y]
                return semilogy(x, y, **kwargs)
            elif isinstance(column_y, (list, tuple)):
                X = list()
                Y = list()
                for i, col_y in enumerate(column_y):
                    assert col_y in self._filer_.columns, f"column_y[{i}] = {col_y} is not valid."
                    X.append(x)
                    Y.append(self._filer_.df[col_y])
                return semilogy(X, Y, num_lines=i+1, **kwargs)
            else:
                raise Exception()
        else:
            raise NotImplementedError()

    def loglog(self, column_x,  column_y, **kwargs):
        """We plot data of 'column_x' (x-axis) and data of 'column_y' (y-axis).

        Parameters
        ----------
        column_x
        column_y
        kwargs

        Returns
        -------

        """
        if RANK != MASTER_RANK: return

        if isinstance(column_x, str):
            assert column_x in self._filer_.columns, f"column_x = {column_x} is not valid."
            x = self._filer_.df[column_x]
            if isinstance(column_y, str):
                # plot one single line.
                assert column_y in self._filer_.columns, f"column_y = {column_y} is not valid."
                y = self._filer_.df[column_y]
                return loglog(x, y, **kwargs)
            elif isinstance(column_y, (list, tuple)):
                X = list()
                Y = list()
                for i, col_y in enumerate(column_y):
                    assert col_y in self._filer_.columns, f"column_y[{i}] = {col_y} is not valid."
                    X.append(x)
                    Y.append(self._filer_.df[col_y])
                return loglog(X, Y, num_lines=i+1, **kwargs)
            else:
                raise Exception()
        else:
            raise NotImplementedError()




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
