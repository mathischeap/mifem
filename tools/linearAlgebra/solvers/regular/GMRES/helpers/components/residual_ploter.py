# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
from components.miscellaneous.timer import MyTimer


def ___gmres_plot_residuals___(residuals, solve_name, scheme_name):
    """

    :param residuals: We plot these residuals.
    :param solve_name: We will save the plot using this name.
    :return:
    """
    assert solve_name is not None, f"to plot residuals, give the solving process a name."

    fig = plt.figure()
    plt.semilogy(residuals)
    plt.title(scheme_name + ': ' + solve_name  + '@' + MyTimer.current_time())
    plt.ylabel('residuals')
    plt.xlabel('iterations')

    if os.path.isdir("__residuals__"):
        pass
    else:
        os.mkdir('__residuals__')

    plt.savefig(f'./__residuals__/'
        f'{scheme_name + "_" + solve_name + "_" + MyTimer.current_time_with_no_special_characters()}.png')

    plt.close(fig)