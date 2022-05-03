# -*- coding: utf-8 -*-


from tools.deprecated.assembler.main import GatheringMatrix

import numpy as np
from root.config.main import rAnk, mAster_rank, cOmm, MPI

from tools.linear_algebra.gathering.regular.matrix.do.main import Gathering_Matrix_DO

class Gathering_Matrix(GatheringMatrix):
    """A gathering matrix is a bunch of gathering vectors.

    Gathering vectors (as values) stored in a dictionary whose keys are mesh element indices.
    So it can not be used as, for example, trace-element-wise gathering matrix.

    :param dict gvd: A dictionary of gathering vectors.
    """
    def __init__(self, gvd, mesh_type=None):
        assert isinstance(gvd, dict), "I need a dictionary of all local gathering vectors."
        self._gvd_ = gvd
        self._GLOBAL_num_dofs_ = None
        self._GLOBAL_len_ = None

        len_i = []
        for i in self:
            gvi = self[i]
            len_i.append(len(gvi))
        len_i = np.array(len_i)
        if np.size(len_i) > 0:
            assert np.all(len_i == len_i[0])
            LOCAL_len = len_i[0]
        else:
            LOCAL_len = -1
        LOCAL_len = cOmm.gather(LOCAL_len, root=mAster_rank)
        if rAnk == mAster_rank:
            shape_1 = np.max(LOCAL_len)
            for i in LOCAL_len:
                if i != shape_1: assert i == -1
        else:
            shape_1 = None
        shape_1 = cOmm.bcast(shape_1, root=mAster_rank)

        self._GLOBAL_shape_ = (self.GLOBAL_len, shape_1)
        self._local_range_ = None
        self._mesh_type_ = mesh_type
        self._do_ = None
        self._freeze_self_()

    def __getitem__(self, i):
        """Return the GV for mesh element #`i`."""
        return self._gvd_[i]

    def __iter__(self):
        """Go through all local mesh elements."""
        for i in self._gvd_:
            yield i

    def __contains__(self, item):
        """If a local mesh element (not a dof) is contained by this GM."""
        return item in self._gvd_

    def __len__(self):
        """local length, equal to the amount of local elements."""
        return len(self._gvd_)

    def __eq__(self, other):
        if self is other:
            RETURN = True
        else:
            if other.__class__.__name__ != 'Gathering_Matrix':
                RETURN = False
            elif len(self) != len(other):
                RETURN = False
            else:
                RETURN = True
                for i in self:
                    if self[i] != other[i]:
                        RETURN = False
                        break

        RETURN = cOmm.allreduce(RETURN, op=MPI.LAND)
        return RETURN

    @property
    def GLOBAL_num_dofs(self):
        """How many dofs in total in all cores."""
        if self._GLOBAL_num_dofs_ is None:
            LOCAL_MAX = [-1,]
            for i in self:
                gv = self[i]
                LOCAL_MAX.append(gv.___PRIVATE_find_max_label___())
            LOCAL_MAX = np.max(LOCAL_MAX)
            self._GLOBAL_num_dofs_ = cOmm.allreduce(LOCAL_MAX, op=MPI.MAX)
            assert self._GLOBAL_num_dofs_ >= 0
            self._GLOBAL_num_dofs_ += 1
        return self._GLOBAL_num_dofs_

    @property
    def GLOBAL_len(self):
        """The global length: the global number of elements."""
        if self._GLOBAL_len_ is None:
            LEN = len(self)
            LEN = cOmm.allreduce(LEN, op=MPI.SUM)
            self._GLOBAL_len_ = LEN
        return self._GLOBAL_len_

    @property
    def GLOBAL_shape(self):
        """Return (global number of elements, local number of dofs)."""
        return self._GLOBAL_shape_

    @property
    def mesh_type(self):
        """"""
        return self._mesh_type_

    @property
    def local_range(self):
        """The range of numbering of dofs in this core. So the dofs of this core are in this range, but
        they do not necessarily fully cover this range.

        If the local dofs are 0,1,2,4,5,12,32,33,35, then the  local_range is (0,36). Note it is not 35.
        So the dofs must be in range(*local_range).
        """
        if self._local_range_ is None:
            if len(self) == 0:
                self._local_range_ = tuple()
            else:
                LOCAL_MIN = list()
                LOCAL_MAX = list()
                for i in self:
                    gv = self[i]
                    LOCAL_MIN.append(gv.___PRIVATE_find_min_label___())
                    LOCAL_MAX.append(gv.___PRIVATE_find_max_label___())
                self._local_range_ = (np.min(LOCAL_MIN), np.max(LOCAL_MAX)+1)
        return self._local_range_


    @property
    def do(self):
        if self._do_ is None:
            self._do_ = Gathering_Matrix_DO(self)
        return self._do_
