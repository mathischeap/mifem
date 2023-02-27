# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/2/2022 4:37 PM
"""
from components.freeze.main import FrozenOnly
from components.assemblers import ___Pr_gathering_matrix_checker___


class VectorDistributor(FrozenOnly):
    """A distributor does the opposite as an assembler."""

    def __init__(self, gathering):
        """"""
        if isinstance(gathering, dict):
            self._indices_ = gathering.keys()
        else:
            self._indices_ = range(len(gathering))
        self._gathering_ = gathering
        self._default_routine_ = 'basic'
        self._num_dofs_ = ___Pr_gathering_matrix_checker___(gathering, self._indices_)
        self._freeze_self_()

    def __call__(self, TBD, routine=None, **kwargs):
        """The output will be a dict or list (depends on gathering) of 1d arrays."""
        if routine is None:
            routine = self._default_routine_
        else:
            pass
        return getattr(self, '___' + routine + "___")(TBD, **kwargs)

    def ___basic___(self, TBD, scheme=1):
        """"""
        if TBD.__class__.__name__ == 'ndarray':
            assert TBD.ndim == 1, f"I need a 1d array."
        else:
            pass
        assert len(TBD) == self._num_dofs_, f"data length wrong."

        if scheme is None:
            scheme = 1
        else:
            pass

        if scheme == 1:

            GM = self._gathering_
            indices = self._indices_

            if isinstance(GM, dict):
                output = dict()
                for i in indices:
                    output[i] = TBD[GM[i]]

            else:
                output = list()
                for i in indices:
                    output.append(TBD[GM[i]])

            return output

        else:
            raise NotImplementedError(f"scheme {scheme} is not implement for "
                                      f"distributor routine: basic.")
