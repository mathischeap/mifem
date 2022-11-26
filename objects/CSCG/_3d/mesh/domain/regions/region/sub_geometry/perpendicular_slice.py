

from components.freeze.main import FrozenOnly

class RegionPerpendicularSlice(FrozenOnly):
    """A perpendicular slice object in the regions.

    The local coordinate of a regions is (r, s, t) = [0,1]^3.

    If r = 0.5, s=[0,1], t=[0, 1] then it means a full slice at (local coordinate) r=0.5.

    """

    def __init__(self, region, r=None, s=None, t=None):
        """"""
        self._region_ = region
        r, s, t, perpendicular_to_axis = self.___PRIVATE_check_x_y_z___(r, s, t)
        self._r_, self._s_, self._t_ = r, s, t
        self._PTA_ = perpendicular_to_axis
        self._freeze_self_()

    @staticmethod
    def ___PRIVATE_check_x_y_z___(r, s, t):
        """"""
        if r is None: r = [0, 1]
        if s is None: s = [0, 1]
        if t is None: t = [0, 1]

        NUM_float, NUM_range = 0, 0
        for i, _ in enumerate([r, s, t]):
            if isinstance(_, (int, float)):
                NUM_float += 1
                assert 0 <= _ <= 1, 'rst'[i] + f'={_} is wrong.'
                perpendicular_to_axis = i
            elif isinstance(_, (tuple, list)):
                NUM_range += 1
                assert len(_) == 2, 'rst'[i] + f'={_} is wrong.'
                S, E = _
                assert 0 <= S <= 1, 'rst'[i] + f'={_} is wrong.'
                assert 0 <= E <= 1, 'rst'[i] + f'={_} is wrong.'
                assert S < E      , 'rst'[i] + f'={_} is wrong.'
            else:
                raise Exception('rst'[i] + f'={_}({_.__class__}) is wrong.')

        assert NUM_float == 1 and NUM_range == 2, f"r, s ,t = {r}, {s}, {t} wrong."
        # noinspection PyUnboundLocalVariable
        perpendicular_to_axis = 'rst'[perpendicular_to_axis]

        return r, s, t, perpendicular_to_axis


    @property
    def r(self):
        """Or x-direction."""
        return self._r_


    @property
    def s(self):
        """Or y-direction."""
        return self._s_

    @property
    def t(self):
        """Or z-direction."""
        return self._t_

    @property
    def perpendicular_to_axis(self):
        return self._PTA_