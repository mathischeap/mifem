# -*- coding: utf-8 -*-


from screws.freeze.main import FrozenOnly
import numpy as np
from itertools import chain



class Gathering_Vector(FrozenOnly):
    """A gathering vector stores the numbering of dofs in one element.

    Indices if the vector represent the local numbering of the dofs.

    :param int i: The element this vector is representing.
    :param gv: The local numbering vector. It must be iterable, like a 1d list, tuple, range, or
        ndarray.
    """
    def __init__(self, i, gv):
        self._i_ = i

        if isinstance(gv, list): # list of integers
            self._gv_ = np.array(gv)

        if isinstance(gv, tuple): # tuple of a series of range object.
            for gvi in gv: assert isinstance(gvi, range), "Tuple of ranges."
            CHAIN = chain(*gv)
            self._gv_ = np.array([_ for _ in CHAIN])
        elif gv.__class__.__name__ == 'ndarray': # 1d numpy.array of integers.
            assert np.ndim(gv) == 1, f"gathering vector needs to be 1d."
            self._gv_ = gv
        elif isinstance(gv, range): # A single range object.
            self._gv_ = gv
        else:
            raise Exception(f'gathering vector type {gv.__class__.__name__} wrong.')
        # do check, we only accept 1d array or range ...
        if self._gv_.__class__.__name__ == 'ndarray':
            assert np.ndim(self._gv_) == 1
        else:
            assert isinstance(self._gv_, range)
        # ...
        self._full_vector_ = None
        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def full_vector(self):
        """We make a full vector for fast usage. self._gv_ can be np.array or a range object."""
        if self._full_vector_ is None:
            if self._gv_.__class__.__name__ == 'ndarray':
                self._full_vector_ = self._gv_
            elif isinstance(self._gv_, range):
                self._full_vector_ = np.array([_ for _ in self])
            else:
                raise Exception()

        return self._full_vector_

    def __getitem__(self, index):
        """Get the ith dof numbering."""
        return self.full_vector[index]

    def __contains__(self, dof):
        """If #i dof is contained by this GV."""
        return dof in self._gv_

    def __iter__(self):
        """Go through all local dofs."""
        for dof in self._gv_:
            yield dof

    def __len__(self):
        """How many local dofs?"""
        return len(self._gv_)

    def __eq__(self, other):
        """If two GVs are equal?"""
        if self is other:
            return True
        else:
            if other.__class__.__name__ != 'Gathering_Vector':
                return False
            elif self.i != other.i:
                return False
            elif len(self) != len(other):
                return False
            else:
                return np.all(self.full_vector == other.full_vector)

    def __repr__(self):
        return 'element#' + str(self.i) + '-' +\
               f"({self.___PRIVATE_find_min_label___()}, " \
               f"{self.___PRIVATE_find_max_label___()})"

    def ___PRIVATE_find_max_label___(self):
        return np.max(self.full_vector)

    def ___PRIVATE_find_min_label___(self):
        return np.min(self.full_vector)

    def index(self, dof):
        """find the index of dof #i. This is like the index function of a list. For example,
            >>> A = [1,2,5,3,4],
            >>> A.index(5)
            2
        """
        WHERE =  np.where(self.full_vector == dof)[0]
        assert len(WHERE) == 1, f"We must only find one dof is numbered i."
        return WHERE[0]