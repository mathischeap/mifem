# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/8 13:05
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class nLS_NonlinearTerms(FrozenOnly):
    """We save the nonlinear terms of a nonlinear system in a list. If the system A is of
    block shape (m, n), we must has m nonlinear terms. Each item of the nonlinear terms must
    be one of (0, None, other).

    Other can be an object of a list or tuple of objects. If itis a list or tuple of objects,
    the sum of the these objects are the nonlinear terms for this block row.

    For example,

    nonlinear terms = [0, 0, 0] : the system is a linear system, actually.
    nonlinear terms = (MDM, 0, 0) : the first equation has an MDM nonlinear term.
    nonlinear terms = (None, [MDM, MDM], 0) : the second equation has two MDM nonlinear terms.

    """
    def __init__(self, nLS):
        """"""
        self._nLS_ = nLS
        self._terms_ = nLS.___nonlinear_terms___

        self._IBD_ = dict()
        self._num_ = 0 # number of nonlinear terms.
        for i, term in enumerate(self._terms_):
            if isinstance(term, (list, tuple)):
                for j, tm in enumerate(term):
                    index = (i, j)
                    assert tm is not None and tm != 0
                    self._IBD_[index] = tm
                    self._num_ += 1
            else:
                if term is None or term == 0:
                    pass
                else:
                    self._num_ += 1
                    self._IBD_[(i,)] = term

        self._freeze_self_()

    def __iter__(self):
        """Go through indices of all nonlinear terms."""
        for index in self._IBD_:
            yield index

    def __getitem__(self, index):
        """Get all nonlinear terms"""
        if isinstance(index, tuple):
            # return a particular nonlinear term
            return self._IBD_[index]
        else:
            # return the nonlinear term(s) for the ith set of equations.
            return self._terms_[index]

    @property
    def num(self):
        """Number of nonlinear terms involved."""
        return self._num_

    def ___PRIVATE_evaluate___(self, unknown_variables):
        """"""
        assert len(unknown_variables) == self._nLS_.shape[1]

        pairs = list()
        for i, uv in enumerate(self._nLS_.unknown_variables):
            pairs.append((uv, unknown_variables[i]))

        nonlinear_f = list()

        for i in range(self._nLS_.shape[0]):
            nfi = list()
            NTs = self[i]
            if NTs == 0 or NTs is None:
                pass
            else:
                if isinstance(NTs, (list, tuple)):
                    for nts in NTs:
                        EWC_v = nts.do.evaluate(pairs)
                        nfi.append(EWC_v)
                else:
                    nfi.append(NTs.do.evaluate(pairs))

            nonlinear_f.append(nfi)

        return nonlinear_f






if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
