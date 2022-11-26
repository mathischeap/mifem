# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/31 12:21 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm

def __matplot_contour_f__(plot_type,
    # data, and linewidth
    x, y,
    # style, color, and labels
    levels=None, num_levels=20, linewidth=1, linestyle=None,
    # config
    usetex=True, saveto=None, corlormap='Dark2',
    # figure
    figsize=(5.5,4), left=0.15, bottom=0.15,
    # title
    title = None, title_size=20, title_pad=12,
    # labels
    xlabel=None, ylabel=None, label_size=16,
    # ticks
    tick_style= 'sci', xticks = None, yticks = None,
    tick_size=16, tick_pad=6, minor_tick_length=4, major_tick_length=8,
    # legend
    legend_size=18, legend_local='best', legend_frame=False,
    ):
    """"""




if __name__ == "__main__":
    # mpiexec -n 4 python components/matplot_wrappers/contour.py
    pass
