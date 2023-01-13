# -*- coding: utf-8 -*-
""" We test the runner with a csv reader.

This is more like an example rather than a test.
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/12/2022 11:21 AM
"""
import os
import sys
if './' not in sys.path:
    sys.path.append('./')

from components.miscellaneous.docstringReaders.numpy_styple import NumpyStyleDocstringReader
import pandas as pd

current_dir = os.path.dirname(__file__)
from tools.run.reader import ParallelMatrix3dInputRunner
from components.miscellaneous.mios import remove, isfile


def csvReader(N, K, c):
    """

    Parameters
    ----------
    N :
    K :
    c :

    Returns
    -------
    a_L2_error :
    b_L2_error :

    """
    if N == 1:
        steps = 20
    else:
        if K == 8:
            steps = 60
        else:
            steps = 40

    filename = current_dir + f"/N{N}K{K}steps{steps}c{c}.csv"

    data = pd.read_csv(filename, index_col=0)
    array = data.to_numpy()
    keys = list(data.keys())

    returns = NumpyStyleDocstringReader(csvReader).Returns

    RT = list()
    for rt in returns:
        RT.append(array[-1, keys.index(rt)])

    return RT


Ns = [[1, 1, 1, ],
      [2, 2, 2, ]]
Ks = [[4, 6, 8, ],
      [4, 6, 8, ]]
cs = [0, ]

pr = ParallelMatrix3dInputRunner(csvReader)

if isfile(current_dir + '/csv_reader.txt'):
    remove(current_dir + '/csv_reader.txt')

pr.iterate(Ns, Ks, cs, writeto=current_dir + '/csv_reader.txt', show_progress=False)

if isfile(current_dir + '/csv_reader.txt'):
    remove(current_dir + '/csv_reader.txt')

try:
    pr.visualize('loglog', 'N', 'a_L2_error', prime='input2', hcp=1, usetex=True,
                 labels=['$N=1$', '$N=2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                 title=False,
                 yticks=[1e-1, 1e-2, 1e-3],
                 xlabel=r'$1/K$',
                 ylabel=r"$\left\| a^h \right\|_{L^2-\mathrm{error}}$",
                 order_text_size=15,
                 plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 1},
                                      1: {'tp': (0.02, 0.2), 'order': 2}},
                 saveto=current_dir + '/a_L2_error.png')

    pr.visualize('loglog', 'N', 'b_L2_error', prime='input2', hcp=1, usetex=True,
                 labels=['$N=1$', '$N=2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                 title=False,
                 yticks=[1e-1, 1e-2, 1e-3, 1e-4],
                 xlabel=r'$1/K$',
                 ylabel=r"$\left\| b^h \right\|_{L^2-\mathrm{error}}$",
                 order_text_size=15,
                 plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 2},
                                      1: {'tp': (0.02, 0.2), 'order': 3}},
                 saveto=current_dir + '/b_L2_error.png')

except RuntimeError:
    pr.visualize('loglog', 'N', 'a_L2_error', prime='input2', hcp=1, usetex=False,
                 labels=['$N=1$', '$N=2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                 title=False,
                 yticks=[1e-1, 1e-2, 1e-3],
                 xlabel=r'$1/K$',
                 order_text_size=15,
                 plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 1},
                                      1: {'tp': (0.02, 0.2), 'order': 2}},
                 saveto=current_dir + '/a_L2_error.png')

    pr.visualize('loglog', 'N', 'b_L2_error', prime='input2', hcp=1, usetex=False,
                 labels=['$N=1$', '$N=2$'],
                 styles=["-s", "-v"],
                 colors=[(0, 0, 0, 1), (0.4, 0.4, 0.4, 1), (0.75, 0.75, 0.75, 1)],
                 title=False,
                 yticks=[1e-1, 1e-2, 1e-3, 1e-4],
                 xlabel=r'$1/K$',
                 order_text_size=15,
                 plot_order_triangle={0: {'tp': (0.02, -0.5), 'order': 2},
                                      1: {'tp': (0.02, 0.2), 'order': 3}},
                 saveto=current_dir + '/b_L2_error.png')


test_csv_reader_passed = 1  # do not comment this.


if __name__ == '__main__':
    # python tests/tools/ParallelMatrix3dInputRunner/csv_reader_test.py
    csvReader(1, 4, 0)
