# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/31 12:41 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from objects.CSCG._2d.discreteFields.base.portion import _2dCSCG_DF_PortionBase
import numpy as np


class _2dCSCG_DF_ScalarPortion(_2dCSCG_DF_PortionBase):
    """A partial of the 2d CSCG region-wise scalar. It again is a region-wise scalar but in each
    region we only save a part of the data."""

    def __init__(self, df):
        """"""
        super(_2dCSCG_DF_ScalarPortion, self).__init__(df)
        self._freeze_self_()

    def xy_range(self, x, y):
        """Get a  portion through x- and y-ranges. For example x = (x0, x1), y=(y0, y1), the
        results will be in rectangle x times y.

        Parameters
        ----------
        x : {list, tuple}
            A list or tuple of two floats.
        y : {list, tuple}
            A list or tuple of two floats.

        Returns
        -------

        """
        assert self._df_.structured, f"xy_range portion only works for structured data."
        assert self._df_.grid is not None, f"xy_range needs grid."

        INDICES = self.___PRIVATE_parse_region_wise_structured_data_indices___(x, y)

        COO = dict()
        VAL = dict()

        for rn in INDICES:
            coo = self._df_.coordinates[rn]
            val = self._df_.values[rn]
            indices = INDICES[rn]
            xInd, yInd = indices
            x0, x1 = xInd
            y0, y1 = yInd

            _coo = list()
            for xy in coo:
                _coo.append(xy[x0:x1, y0:y1])
            _val = list()
            for v in val:
                _val.append(v[x0:x1, y0:y1])

            COO[rn] = _coo
            VAL[rn] = _val

        CLASS = self._df_.__class__

        range_df = CLASS(self._df_.mesh, COO, VAL, 'xy-range-of-' + self._df_.name,
                         structured=False, grid=None)

        return range_df


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/discreteFields/scalar/portion.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0, bounds=([-1, 1], [-1, 1]))([5, 5])
    # mesh = MeshGenerator('chp1',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)
    ES = ExactSolutionSelector(mesh)('sL:sincos1')
    f0 = FC('0-f-o', hybrid=False)
    f0.CF = ES.potential

    f0.CF.current_time = 0
    f0.discretize()

    x = np.linspace(-1, 1, 200)

    ds = f0.reconstruct.discrete_scalar([x, x])

    portion = ds.portion

    pdf = portion.xy_range([-0.25, 0.5], [0.75, 2])

    pdf.visualize.matplot.contour()
