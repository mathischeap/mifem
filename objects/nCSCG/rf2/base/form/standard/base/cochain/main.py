# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 12:24 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from screws.exceptions import _2nCSCG_RF2_SignatureDisMatchError


class nCSCG_RF2_StandardFormCochainBase(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._local_ = None
        self.___Pr_reset_cache___()
        self._freeze_self_()

    def ___Pr_reset_cache___(self): # when we have new representations of cochain, add their cache here
        """"""

    def ___Pr_RGW_cochain___(self):
        """Must have this for saving a form."""
        raise NotImplementedError()

    @property
    def local(self):
        if self._local_ is None:
            raise Exception(f"local cochain is None.")

        if self._local_.signature == self._f_.signature:
            return self._local_
        else:
            raise _2nCSCG_RF2_SignatureDisMatchError('Signature dis-match, renew cochain first.')

    @local.setter
    def local(self, local):
        """"""
        assert local.___Pr_is_nCSCG_RF2_Standard_LC___
        assert local._signature_ == self._f_.signature, f"signature does not match." # update the signature.
        self._local_ = local








if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
