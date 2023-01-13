# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/31 12:21 PM
"""
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import numpy as np
from root.config.main import RANK, SECRETARY_RANK, COMM


def ___set_contour_levels___(v, num_levels):
    """ We set the `num_levels` according to the values `v` to be plotted. """
    MAX = list()
    MIN = list()
    for rn in v:
        v_ = v[rn][0]
        MAX.append(np.max(v_))
        MIN.append(np.min(v_))

    MINv = np.min(MIN)
    MAXv = np.max(MAX)
    if MINv == MAXv:
        if MINv == 0:
            MAXv = 0.1
        else:
            if MINv > 0:
                MAXv = MINv * 1.01
            else:
                MAXv = MINv * 0.99
        num_levels = 2
    levels = np.linspace(MINv, MAXv, num_levels)
    return levels


def contour(
        xy, v,
        levels=None, num_levels=20, linewidth=1, linestyles=None,
        usetex=False,
        figsize=(6, 5),
        colormap='coolwarm', color=None,  # color will override color map, only works in contour.
        axis_on=True, tick_on=True,
        xlim=None, ylim=None,
        show_colorbar=True,
        colorbar_label=None, colorbar_orientation='vertical', colorbar_aspect=20,
        colorbar_labelsize=13, colorbar_extend='both',
        xticks=None, yticks=None,
        tick_size=15, tick_pad=6, minor_tick_length=4, major_tick_length=8,
        label_size=18,
        title=None,
        saveto=None, dpi=210,
        plot_type='contour',
        colorbar_only=False,
        negative_linestyles=None,
):
    """

    Parameters
    ----------
    xy: {dict, 3d-array}
    v: {dict, 3d-array}
    levels
    num_levels :
        Only applies when `levels` is None.

    linewidth : {float, int, array-like}
        Only applies to contour.

        The line width of the contour lines.

        If a number, all levels will be plotted with this linewidth.

        If a sequence, the levels in ascending order will be plotted with the linewidths in
        the order specified.

        If None, this falls back to rcParams["lines.linewidth"] (default: 1).

    linestyles
        Only applies to contour.

        If linestyles is None, the default is 'solid' unless the lines are monochrome.
        In that case, negative contours will take their linestyles from
        rcParams["contour.negative_linestyle"] (default: 'dashed') setting.

        linestyles can also be an iterable of the above strings specifying a
        set of linestyles to be used. If this iterable is shorter than the number of
        contour levels it will be repeated as necessary.

    usetex : bool, optional, default: False
        If we use LaTeX render?

    colormap
    color : (None, str), optional
        Only applies to contour.

        When  it is not None, it overrides the colormap.

    axis_on :
        Show frame.
    tick_on :
        Show ticks.

        Only applies when `axis_on` is True.

    xlim
    ylim
    show_colorbar
    figsize
    colorbar_label
    colorbar_orientation
    colorbar_aspect
    colorbar_labelsize
    colorbar_extend
    xticks
    yticks
    tick_size
    tick_pad
    minor_tick_length
    major_tick_length
    label_size
    title
    saveto
    dpi
    plot_type : str, optional
        'contourf' or 'contour'.
    colorbar_only: bool
    negative_linestyles: {None, 'solid', 'dashed', 'dashdot', 'dotted'}, optional

    Returns
    -------

    """
    # ... parse data
    if isinstance(xy, dict):
        assert isinstance(v, dict), f"xy is dict, value must be dict as well."

    elif isinstance(xy, tuple) and len(xy) == 2:
        x, y = xy
        if x.__class__.__name__ == 'ndarray':
            assert y.__class__.__name__ == 'ndarray', f"x is array, y must be a array as well."
            assert v.__class__.__name__ == 'ndarray', f"x and y are arrays, value must be a array as well."

            assert np.shape(x) == np.shape(y) == np.shape(v), f"x, y, v shape dis-match."

            XY = dict()
            V = dict()
            if np.ndim(x) == 2:
                XY[f'{RANK}-block0'] = [x, y]
                V[f'{RANK}-block0'] = v
            elif np.ndim(x) == 3:
                for i, xi in enumerate(x):
                    XY[f'{RANK}-block{i}'] = [xi, y[i]]
                    V[f'{RANK}-block{i}'] = v[i]
            else:
                raise Exception(f"can only plot 2d or 3d data.")

            xy = XY
            v = V
    else:
        raise NotImplementedError(f"cannot handle xy ({xy.__class__.__name__}).")

    # ... check data
    assert isinstance(xy, dict) and isinstance(v, dict), f"we should have put all data into dict already."
    for block_name in v:
        vbn = v[block_name]
        if isinstance(vbn, (list, tuple)) and len(vbn) == 1:
            v[block_name] = vbn[0]

    # ... gather data

    xy = COMM.gather(xy, root=SECRETARY_RANK)
    v = COMM.gather(v, root=SECRETARY_RANK)

    if RANK != SECRETARY_RANK:
        return

    _ = dict()
    for __ in xy:
        _.update(__)
    xy = _

    _ = dict()
    for __ in v:
        _.update(__)
    v = _

    # ----------- prepare levels -------------------------------------------------------------1
    if levels is None:
        levels = ___set_contour_levels___(v, num_levels)
    else:
        pass

    # ------ config plt ----------------------------------------------------------------------1
    if saveto is not None:
        matplotlib.use('Agg')

    plt.rcParams.update({
        "text.usetex": usetex,
        "font.family": "DejaVu sans",
        # "font.serif": "Times New Roman",
    })

    plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
    if colormap is not None:
        plt.rcParams['image.cmap'] = colormap
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')

    # ---------- set labels, ticks, frame and so on ------------------------------------------1
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if axis_on:
        if tick_on:
            if xticks is not None:
                plt.xticks(xticks)
            if yticks is not None:
                plt.yticks(yticks)
            plt.tick_params(which='both', labeltop=False, labelright=False, top=False, right=False)
            plt.tick_params(axis='both', which='minor', direction='out', length=minor_tick_length)
            plt.tick_params(axis='both', which='major', direction='out', length=major_tick_length)
            plt.tick_params(axis='both', which='both', labelsize=tick_size)
            plt.tick_params(axis='x', which='both', pad=tick_pad)
            plt.tick_params(axis='y', which='both', pad=tick_pad)
            plt.xlabel('$x$', fontsize=label_size)
            plt.ylabel('$y$', fontsize=label_size)
        else:
            plt.tick_params(
                which='both',
                labeltop=False, labelright=False, labelleft=False, labelbottom=False,
                top=False, right=False, left=False, bottom=False
            )
    else:
        ax.set_axis_off()

    # ------------- contour or contourf plot -------------------------------------------------1
    for block_name in v:
        if plot_type == 'contour':
            if color is None:
                plt.contour(*xy[block_name], v[block_name],
                            levels=levels, linewidths=linewidth, linestyles=linestyles)
            else:
                plt.contour(*xy[block_name], v[block_name], negative_linestyles=negative_linestyles,
                            colors=color, levels=levels, linewidths=linewidth, linestyles=linestyles)

        elif plot_type == 'contourf':
            VAL = v[block_name]
            VAL[VAL > levels[-1]] = levels[-1]
            VAL[VAL < levels[0]] = levels[0]
            plt.contourf(*xy[block_name], VAL, levels=levels)
        else:
            raise Exception(f"plot_type={plot_type} is wrong. Should be one of ('contour', 'contourf')")

    # --------------- title ------------------------------------------------------------------1
    if title is None:
        pass  # default title: No title.
    elif title is False:
        pass
    else:
        plt.title(title)

    # ------------------------------ color bar ----------------------------------------------1
    if (plot_type == 'contour' and color is None and show_colorbar) or (plot_type == 'contourf' and show_colorbar):
        mappable = cm.ScalarMappable()
        mappable.set_array(np.array(levels))
        cb = plt.colorbar(mappable, ax=ax,
                          extend=colorbar_extend,
                          aspect=colorbar_aspect,
                          orientation=colorbar_orientation)

        if colorbar_label is not None:
            cb.set_label(colorbar_label, labelpad=10, size=15)

        cb.ax.tick_params(labelsize=colorbar_labelsize)

    else:
        pass

    # --------------------- save to ----------------------------------------------------------1
    if colorbar_only:
        ax.remove()
    else:
        pass

    if saveto is None:
        plt.show()
    else:
        if saveto[-4:] == '.pdf':
            plt.savefig(saveto, bbox_inches='tight')
        else:
            plt.savefig(saveto, dpi=dpi, bbox_inches='tight')
    plt.close()

    # ========================================================================================1
    return fig


def contourf(*are, **kwargs):
    """"""
    return contour(*are, **kwargs, plot_type='contourf')
