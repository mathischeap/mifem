
from screws.freeze.main import FrozenOnly
import numpy as np

class _3dCSCG_ADF_SF_Error(FrozenOnly):
    """This is one of the key properties of a algebraic dual form.

    To perform the coboundary, now we need the help of another dual trace form.
    """
    def __init__(self, dsf):
        self._dsf_ = dsf
        self._freeze_self_()


    def L(self, *args, **kwargs):
        """"""
        return self._dsf_.prime.error.L(*args, **kwargs)


    def dH(self, dt, dfunc, time, n=1):
        """

        Parameters
        ----------
        n : int
            We compute the `dual{H}^n`-error.
        dt :
            The dual trace form.

        Returns
        -------

        """

        dual_d_form = self._dsf_.coboundary(dt)

        dual_d_form.prime.TW.func.body = dfunc
        dual_d_form.prime.TW.do.push_func_to_instant(time)

        self_error_L2 = self.L()
        ddf_error_L2 = dual_d_form.error.L()

        if n == 1:
            return np.sqrt(self_error_L2**2 + ddf_error_L2**2)
        else:
            raise NotImplementedError()