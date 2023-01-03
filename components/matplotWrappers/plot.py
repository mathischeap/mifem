# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/16 5:12 PM
"""
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm


def __matplot__(
    plot_type,
    # data, and linewidth
    x, y, num_lines=1, linewidth=1.2,
    # style, color, and labels
    style=None, color=None, label=False,
    styles=None, colors=None, labels=None,
    # config
    usetex=True, saveto=None, corlormap='Dark2',
    # figure
    figsize=(5.5, 4), left=0.15, bottom=0.15,
    # title
    title=None, title_size=20, title_pad=12,
    # labels
    xlabel=None, ylabel=None, label_size=16,
    # ticks
    tick_style='sci', xticks=None, yticks=None,
    tick_size=16, tick_pad=6, minor_tick_length=4, major_tick_length=8,
    # legend
    legend_size=18, legend_local='best', legend_frame=False,
):
    """

    Parameters
    ----------
    plot_type
    x :
        If num_lines > 1, we plot (x[i], y[i]) for i in range(num_lines).
    y :
        If num_lines > 1, we plot (x[i], y[i]) for i in range(num_lines).
    num_lines
    style
    color
    label : bool, None, str
        The arg affects when `num_lines` == 1.

        If it is False (default): we will turn off the label.
        If it is None: we will use a default label.
        If it is a str: we will it as the label of the single line.

    styles :
    colors :
    labels : list, tuple, bool
    linewidth
    usetex
    saveto
    corlormap
    figsize
    left
    bottom
    title
    title_size
    title_pad
    xlabel
    ylabel
    label_size
    tick_style : str
        {'sci', 'scientific', 'plain'}
    xticks
    yticks
    tick_size : int
        The font size of the ticks.
    tick_pad : int
        The gap size between the ticks and the axes.
    minor_tick_length
    major_tick_length
    legend_size
    legend_local
    legend_frame

    Returns
    -------

    """
    # - config matplotlib -------------------------------------------------------------------------1
    if saveto is not None:
        matplotlib.use('Agg')
    plt.rc('text', usetex=usetex)
    if usetex:
        plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
    _, ax = plt.subplots(figsize=figsize)
    plt.gcf().subplots_adjust(left=left)
    plt.gcf().subplots_adjust(bottom=bottom)

    # --------- default styles, colors, and labels -------------------------------------------------

    if num_lines == 1:  # single line
        if style is None:
            style = '-^'
        if color is None:
            color = 'darkgray'
        assert label in (None, False) or isinstance(label, str), f"label must be a str, False, or None."
        if label is None:
            label = r'$\mathrm{line}\#0$'

    elif num_lines > 1:  # multiple lines
        if styles is None:
            styles = ('-^', '-x', '-o', '-s', '-v', '-*', '-8', '->', '-p', '-H', '-h', '-D', '-d', '-P') * 5

        if colors is None:
            color = cm.get_cmap(corlormap, num_lines)
            colors = []
            for j in range(num_lines):
                colors.append(color(j))

        if labels is None:
            labels = [r'$\mathrm{line}\#' + str(_) + '$' for _ in range(num_lines)]

        assert isinstance(styles, (list, tuple)) and len(styles) == num_lines, \
            f"put correct amount of styles in list pls."
        assert isinstance(colors, (list, tuple)) and len(styles) == num_lines, \
            f"put correct amount of colors in list pls."

        if labels is False:
            pass
        else:
            assert isinstance(labels, (list, tuple)), f"put labels in list pls."
            assert len(labels) == num_lines, f"I need {num_lines} labels, now I get {len(labels)}."
            for i, lab in enumerate(labels):
                assert isinstance(lab, str), f"labels[{i}] = {lab} is not str."

    else:
        raise Exception()

    # - do the plot -------------------------------------------------------------------------------1
    plotter = getattr(plt, plot_type)

    if num_lines == 1:
        if label is not False:
            plotter(x, y, style, color=color, label=label, linewidth=linewidth)
        else:
            plotter(x, y, style, color=color, linewidth=linewidth)

    else:

        if labels is False:
            for i in range(num_lines):
                plotter(x[i], y[i], styles[i], color=colors[i], linewidth=linewidth)
        else:
            for i in range(num_lines):
                plotter(x[i], y[i], styles[i], color=colors[i], label=labels[i], linewidth=linewidth)

    # ------ customize figure ---------------------------------------------------------------------1
    if xticks is not None:
        plt.xticks(xticks)
    if yticks is not None:
        plt.yticks(yticks)
    plt.tick_params(which='both', labeltop=False, labelright=False, top=True, right=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=minor_tick_length)
    plt.tick_params(axis='both', which='major', direction='in', length=major_tick_length)
    plt.tick_params(axis='both', which='both', labelsize=tick_size)
    plt.tick_params(axis='x', which='both', pad=tick_pad)
    plt.tick_params(axis='y', which='both', pad=tick_pad)

    if plot_type in ('semilogy', 'loglog'):
        pass
    else:
        plt.ticklabel_format(style=tick_style, axis='y', scilimits=(0, 0))
        tx = ax.yaxis.get_offset_text()
        tx.set_fontsize(tick_size)

    if xlabel is not None:
        plt.xlabel(xlabel, fontsize=label_size)
    if ylabel is not None:
        plt.ylabel(ylabel, fontsize=label_size)
    if title is not None:
        plt.title(r'' + title, fontsize=title_size, pad=title_pad)

    # ----legend ----------------------------------------------------------------------------------1
    if num_lines == 1 and label is False:
        # turn off legend when num_line==1 and label is False
        pass
    elif num_lines > 1 and labels is False:
        # turn off legend when num_line > 1 and labels is False
        pass
    else:
        plt.legend(fontsize=legend_size, loc=legend_local, frameon=legend_frame)

    # ---------------- save the figure ------------------------------------------------------------1
    plt.tight_layout()
    if saveto is not None and saveto != '':
        plt.savefig(saveto, bbox_inches='tight')
    else:
        plt.show()
    plt.close()
    # =============================================================================================1
    return


def plot(*args, **kwargs):
    return __matplot__('plot', *args, **kwargs)


def semilogy(*args, **kwargs):
    return __matplot__('semilogy', *args, **kwargs)


def loglog(*args, **kwargs):
    return __matplot__('loglog', *args, **kwargs)
