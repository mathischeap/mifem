from screws.freeze.base import FrozenOnly
from objects.CSCG.base.forms.base.time_wise.BC.do import CSCG_FTWBC_DO



class CSCG_Form_TimeWise_BC(FrozenOnly):
    """Note that the TW.BC of a form will not be saved!!!!.
    """
    def __init__(self, tw):
        self._tw_ = tw
        self._f_ = tw._f_
        self._ndim_ = self._f_.ndim
        self._DO_ = CSCG_FTWBC_DO(self)
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        """Clear cache."""
        self._body_ = None
        self._ES_ = None
        self._ES_variable_name_ = None

    def ___DO_set_BC_body_as___(self, scalar_vector):
        """
        We use this function such that we can set BC.

        :param scalar_vector: scalar, vector
        :return:
        """
        self.body = scalar_vector  # has to pass the checker.


    @property
    def body(self):
        """The function body."""
        return self._body_

    @body.setter
    def body(self, body):
        self._f_.___PRIVATE_TW_BC_body_checker___(body)  # has to pass the checker.
        self._body_ = body
        self._ES_ = None
        self._ES_variable_name_ = None

    @property
    def ftype(self):
        """(str) The function type."""
        return self.body.ftype

    @property
    def DO(self):
        return self._DO_
