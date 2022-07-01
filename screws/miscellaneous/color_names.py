# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/21 3:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import matplotlib.colors as mcolors


def css_color_names():
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(color))), name)
                    for name, color in mcolors.CSS4_COLORS.items())
    names = [name for hsv, name in by_hsv]
    return names

def base_color_names():
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(color))), name)
                    for name, color in mcolors.BASE_COLORS.items())
    names = [name for hsv, name in by_hsv]
    return names

def tableau_color_names():
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(color))), name)
                    for name, color in mcolors.TABLEAU_COLORS.items())
    names = [name for hsv, name in by_hsv]
    return names