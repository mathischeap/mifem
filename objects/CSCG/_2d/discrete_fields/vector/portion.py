# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/31 6:12 PM
"""
import sys
import numpy as np

if './' not in sys.path: sys.path.append('./')
from objects.CSCG._2d.discrete_fields.base.portion import _2dCSCG_DF_PortionBase


class _2dCSCG_DF_VectorPortion(_2dCSCG_DF_PortionBase):
    """A partial of the 2d CSCG region-wise scalar. It again is a region-wise scalar but in each
    region we only save a part of the data."""

    def __init__(self, df):
        """"""
        super(_2dCSCG_DF_VectorPortion, self).__init__(df)
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
    # mpiexec -n 4 python objects/CSCG/_2d/discrete_fields/vector/portion.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    # mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    mesh = MeshGenerator('chp2')([10,10])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    u = FC('1-f-o', is_hybrid=True)

    u.TW.func.do.set_func_body_as(ES, 'velocity')
    u.TW.current_time = 0
    u.TW.do.push_all_to_instant()
    u.discretize()

    r = np.linspace(-1,1,20)
    s = np.linspace(-1,1,20)

    dv = u.reconstruct.discrete_vector([r, s])
    x = [-1., -0]
    y = [-1., -0]

    portion = dv.portion.xy_range(x, y)
    portion.visualize.matplot.quiver(scale=100)
