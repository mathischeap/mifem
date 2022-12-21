# -*- coding: utf-8 -*-
""""""
from components.freeze.base import FrozenOnly
from root.config.main import COMM, MPI, np
from tools.elementwiseCache.gathering.irregular.ir_matrix.do.main import iR_Gathering_Matrix_DO


class iR_Gathering_Matrix(FrozenOnly):
    """Irregular gathering matrix: in each row, the local dofs could be different!"""
    def __init__(self, gvd, mesh_type=None):
        assert isinstance(gvd, dict), f"must be a dict of Gathering_vectors."
        self._gvd_ = gvd
        self._GLOBAL_num_dofs_ = None
        self._GLOBAL_len_ = None

        for key in gvd:
            assert key == gvd[key].i

        self._mesh_type_ = mesh_type
        self._local_range_ = None
        self._do_ = iR_Gathering_Matrix_DO(self)
        self._stamp_ = str(id(self))
        self._freeze_self_()

    @property
    def stamp(self):
        return self._stamp_

    @property
    def ___Pr_IS_regular___(self):
        return False

    def __repr__(self):
        """"""
        return "irGM-"


    def __getitem__(self, rp):
        return self._gvd_[rp]

    def __iter__(self):
        """Go through all local mesh elements."""
        for rp in self._gvd_:
            yield rp

    def __contains__(self, rp):
        """If a local mesh element (not a dof) is contained by this GM."""
        return rp in self._gvd_

    def __len__(self):
        """local length, equal to the amount of local elements."""
        return len(self._gvd_)

    def __eq__(self, other):
        if self is other:
            RETURN = True
        else:
            if other.__class__.__name__ != 'iR_Gathering_Matrix':
                RETURN = False
            elif len(self) != len(other):
                RETURN = False
            else:
                RETURN = True
                for rp in self:
                    if self[rp] != other[rp]:
                        RETURN = False
                        break

        RETURN = COMM.allreduce(RETURN, op=MPI.LAND)
        return RETURN

    @property
    def global_num_dofs(self):
        """How many dofs in total in all cores."""
        if self._GLOBAL_num_dofs_ is None:
            LOCAL_MAX = [-1, ]
            for i in self:
                gv = self[i]
                LOCAL_MAX.append(gv.___PRIVATE_find_max_label___())
            LOCAL_MAX = np.max(LOCAL_MAX)
            self._GLOBAL_num_dofs_ = COMM.allreduce(LOCAL_MAX, op=MPI.MAX)
            assert self._GLOBAL_num_dofs_ >= 0
            self._GLOBAL_num_dofs_ += 1
        return self._GLOBAL_num_dofs_

    @property
    def global_len(self):
        """The global length: the global number of elements."""
        if self._GLOBAL_len_ is None:
            LEN = len(self)
            LEN = COMM.allreduce(LEN, op=MPI.SUM)
            self._GLOBAL_len_ = LEN
        return self._GLOBAL_len_

    @property
    def global_shape(self):
        """Return (global number of elements, local number of dofs)."""
        raise Exception(f"irregular GM has no GLOBAL shape.")


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
        return self._do_
