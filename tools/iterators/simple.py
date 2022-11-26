# -*- coding: utf-8 -*-
"""Time-wise iterator.
"""

from root.config.main import *
from tools.iterators.base.main import Iterator



class SimpleIterator(Iterator):
    """It is 'simple' because we take `dt` is constant!

    :param t0:
    :param dt:
    :param max_steps:
    :param max_time:
    :param kwargs: These kwargs are those passed to the parent: `Iterator`.
    """
    def __init__(self, t0=0, dt=None, max_steps=None, max_time=None, **kwargs):
        CHECK_SAME_IN_ALL_CORES(t0, dt, max_steps, max_time, kwargs)
        self.___PRIVATE_parse_time_max_steps_max_time___(t0, dt, max_steps, max_time)
        super().__init__(**kwargs)

    def __next__(self):
        outputs = self._solver_(self.t, self.t+self.dt)
        # SimpleIterator use constant dt. Updating dt is necessary for all particular iterators.
        self.dt = self.___PRIVATE_update_dt___() # this dt will be used for the next call of the solver
        return outputs

    def ___PRIVATE_update_dt___(self):
        """
        For SimpleIterator, dt is always the same. But we have this method to remind us that we need
        the update method for all iterators.

        :return:
        """
        return self.dt

    def ___PRIVATE_parse_time_max_steps_max_time___(self, t0, dt, max_steps, max_time):
        """ t0: must given! dt > max_steps > max_time."""
        self._t0_ = t0

        if dt is not None:
            assert dt > 0, " <Simple> : dt={} wrong, should > 0.".format(dt)
            self.dt = dt
            if max_steps is not None:
                assert isinstance(max_steps, int) and max_steps > 0, \
                    " <Simple> : max_steps={} wrong.".format(max_steps)
                self._max_steps_ = max_steps
                self._max_time_ = self.t0 + self.dt * self.max_steps
            else:
                if max_time is not None:
                    assert isinstance(max_time, (int, float)) and max_time > self.t0, \
                        " <Simple> : max_time={} wrong!".format(max_time)
                    self._max_steps_ = int(np.ceil((max_time - self.t0) / self.dt))
                else:
                    self._max_steps_ = 9999999
                self._max_time_ = self.t0 + self.dt * self.max_steps

        else:
            assert max_steps is not None, " <Simple> : when dt=None, max_steps cannot be None!"
            assert max_time is not None, " <Simple> : when dt=None, max_time cannot be None!"
            assert isinstance(max_time, (int, float)) and max_time > self.t0, \
                " <Simple> : max_time={} wrong!".format(max_time)
            assert isinstance(max_steps, int) and max_steps > 0, \
                " <Simple> : max_steps={} wrong.".format(max_steps)
            self._max_steps_ = max_steps
            self._max_time_ = max_time
            self.dt = (self.max_time - self.t0) / self.max_steps

        assert self.dt >= 0, " <Simple> : dt wrong!"
        assert isinstance(self.max_steps, int) and self.max_steps > 0, " <Simple> "
        assert self.t0 + self.dt * self.max_steps == self.max_time, " <Simple> : max_time wrong!"