

"""
A perpendicular slice object in the region.

"""

from SCREWS.frozen import FrozenOnly


class RegionPerpendicularSlice(FrozenOnly):
    """
    The local coordinate of a region is (r, t, s) = [0,1]^3.

    If r = 0.5, t=[0,1], s=[0, 1] then it means a full slice at (local coordinate) t=0.5.

    """

    def __init__(self, region, r=None, t=None, s=None):
        """"""
        self._region_ = region
        r, t, s, perpendicular_to_axis = self.___PRIVATE_check_x_y_z___(r, t, s)
        self._r_, self._t_, self._s_ = r, t, s
        self._PTA_ = perpendicular_to_axis
        self._freeze_self_()

    @staticmethod
    def ___PRIVATE_check_x_y_z___(r, t, s):
        """"""
        if r is None: r = [0, 1]
        if t is None: t = [0, 1]
        if s is None: s = [0, 1]

        NUM_float, NUM_range = 0, 0
        for i, _ in enumerate([r, t, s]):
            if isinstance(_, (int, float)):
                NUM_float += 1
                assert 0 <= _ <= 1, 'rts'[i] + f'={_} is wrong.'
                perpendicular_to_axis = i
            elif isinstance(_, (tuple, list)):
                NUM_range += 1
                assert len(_) == 2, 'rts'[i] + f'={_} is wrong.'
                S, E = _
                assert 0 <= S <= 1, 'rts'[i] + f'={_} is wrong.'
                assert 0 <= E <= 1, 'rts'[i] + f'={_} is wrong.'
                assert S < E, 'rts'[i] + f'={_} is wrong.'
            else:
                raise Exception('rts'[i] + f'={_} is wrong.')

        assert NUM_float == 1 and NUM_range == 2, f"r, t,s = {r}, {t}, {s} wrong."
        # noinspection PyUnboundLocalVariable
        perpendicular_to_axis = 'rts'[perpendicular_to_axis]

        return r, t, s, perpendicular_to_axis


    # properties ----------------------- BELOW -------------------------------------------------

    @property
    def r(self):
        return self._r_

    @property
    def t(self):
        return self._t_

    @property
    def s(self):
        return self._s_

    @property
    def perpendicular_to_axis(self):
        return self._PTA_