# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/31 3:24 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2dCSCG_DF_PortionBase(FrozenOnly):
    """"""

    def __init__(self, df):
        """"""
        self._df_ =  df
        self._freeze_self_()


    def ___PRIVATE_parse_region_wise_structured_data_indices___(self, x, y):
        x0, x1 = x
        y0, y1 = y
        assert x0 < x1,  f"x={x} wrong, x0 must be lower than x1."
        assert y0 < y1,  f"y={y} wrong, y0 must be lower than y1."

        assert self._df_.structured, f"xy_range portion only works for structured data."
        assert self._df_.grid is not None, f"xy_range needs grid."
        mesh = self._df_.mesh

        INDICES = dict()

        for rn in self._df_.regions:
            x0, x1 = x # renew initial x0, x1 for each region
            y0, y1 = y # renew initial y0, y1 for each region

            region = mesh.domain.regions[rn]
            rItp = region.interpolation

            x_lim = rItp.mapping_Xr_at_s0([0, 1])
            y_lim = rItp.mapping_Ys_at_r0([0, 1])

            if x0 >= x_lim[1] or x1 <= x_lim[0] or y0 >= y_lim[1] or y1 <= y_lim[0]:
                # this region has no business with the portion data.
                pass

            else:
                if x0 < x_lim[0]: x0 = x_lim[0]
                if x1 > x_lim[1]: x1 = x_lim[1]
                if y0 < y_lim[0]: y0 = y_lim[0]
                if y1 > y_lim[1]: y1 = y_lim[1]
                assert x0 < x1,  f"trivial check, must be"
                assert y0 < y1,  f"trivial check, must be"

                r, s = self._df_.grid[rn]
                R = rItp.mapping_Xr_at_s0(r)
                S = rItp.mapping_Ys_at_r0(s)

                if x0 > R[-1] or x1 < R[0] or y0 > S[-1] or y1 < S[0]:
                    # although x0, x1, y0, y1 are valid, but no data in this range
                    pass
                else:

                    x_ind_l = None
                    for i, r in enumerate(R):
                        if r >= x0:
                            x_ind_l = i
                            break

                    x_ind_u = None
                    for i, r in enumerate(R):
                        if r > x1:
                            x_ind_u = i
                            break
                    if x_ind_u is None:
                        x_ind_u = len(R)

                    y_ind_l = None
                    for j, s in enumerate(S):
                        if s >= y0:
                            y_ind_l = j
                            break

                    y_ind_u = None
                    for j, s in enumerate(S):
                        if s > y1:
                            y_ind_u = j
                            break
                    if y_ind_u is None:
                        y_ind_u = len(S)

                    if x_ind_l == x_ind_u or y_ind_l == y_ind_u:
                        # the grid is too coarse and have no data in the range.
                        pass
                    else:
                        assert x_ind_l < x_ind_u and y_ind_l < y_ind_u, f"trivial check."

                        xInd_bounds = [x_ind_l, x_ind_u]
                        yInd_bounds = [y_ind_l, y_ind_u]

                        INDICES[rn] = (xInd_bounds, yInd_bounds)

        return INDICES

