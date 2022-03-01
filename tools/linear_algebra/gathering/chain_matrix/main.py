# -*- coding: utf-8 -*-
""""""

import sys
if './' not in sys.path: sys.path.append('../')

from screws.frozen import FrozenOnly
from root.config import *
from tools.linear_algebra.gathering.chain_matrix.find import ___Chain_Gathering_Matrix_FIND___

from tools.linear_algebra.gathering.vector import Gathering_Vector
from tools.linear_algebra.gathering.matrix import Gathering_Matrix


class Chain_Gathering_Matrix(FrozenOnly):
    """We chain some Gathering_Matrix together, for example, chain (GM0, GM1, GM2). Then in GM1, the numbering does not
    start with 0 anymore, it starts with ``GM0.GLOBAL_num_dofs``. As for GM2, its numbering now starts with
    ``GM0.GLOBAL_num_dofs + GM1.GLOBAL_num_dofs``.

    The idea is not just chain them, we also provide a mechanism to save memory, we do not make a new Gathering_Matrix
    for this chained Chain_Gathering_Matrix.
    """
    def __init__(self, GMs):
        """
        GMs can be a gathering matrix or a chain gathering matrix. Or a list of them.

        :param GMs:
        """

        if GMs.__class__.__name__ in ('Gathering_Matrix', 'Chain_Gathering_Matrix,'):
            GMs = [GMs,]

        mesh_types = list()

        assert isinstance(GMs, (list, tuple)) and len(GMs) >= 1, \
            "we need a list or tuple of at least one gathering matrices."

        NEW_GMs = list()
        for i, gm in enumerate(GMs):

            mesh_types.append(gm.mesh_type)
            assert mesh_types[i] == mesh_types[0], f"mesh_types[{i}]: {mesh_types[i]} is different from mesh_types[0]: {mesh_types[0]}"


            if gm.__class__.__name__ == 'Gathering_Matrix':
                NEW_GMs.append(gm)
            elif gm.__class__.__name__ == 'Chain_Gathering_Matrix':
                NEW_GMs.extend(gm.GMs)
            else:
                raise Exception(f'GMs[{i}] is a {gm.__class__.__name__}, '
                                f'not a Gathering_Matrix or Chain_Gathering_Matrix.')


        self._mesh_type_ = mesh_types[0]


        GMs = NEW_GMs

        self.___NUM___ = len(GMs)
        # ... check GM in GMs are actually representing same local elements.
        if self.___NUM___ == 1:
            pass
        else:
            gm0 = GMs[0]
            LEN = len(gm0)
            for _, gm in enumerate(GMs[1:]):
                assert len(gm) == LEN, f"GMs[{_+1}] has different length comparing to GMs[0]."
                for i in gm0:
                    assert i in gm, f"GMs[{_+1}] represent different elements comparing to GMs[0]."

        To_Be_Added = [0,]
        for gm in GMs[:-1]:
            GLOBAL_num_dofs = gm.GLOBAL_num_dofs
            # noinspection PyUnresolvedReferences
            To_Be_Added.append(To_Be_Added[-1]+GLOBAL_num_dofs)

        self._To_Be_Added_ = To_Be_Added
        self._GLOBAL_num_dofs_ = To_Be_Added[-1] + GMs[-1].GLOBAL_num_dofs
        self._GMs_ = GMs
        self._NUM_GMs_ = len(GMs)
        self._dofs_distribution_ = [_.GLOBAL_num_dofs for _ in self.GMs]
        self._local_dofs_distribution_ = [_.GLOBAL_shape[1] for _ in self.GMs]
        self._local_ranges_ = None
        self._FIND_ = ___Chain_Gathering_Matrix_FIND___(self)
        if len(self) > 0:
            for i in self:
                assert sum(self.GV_LENs) == len(self[i]), "A safety check."
                break
            self._LOCAL_To_Be_Added_ = [0,] # when len(self) == 0, we do not even initialize self._LOCAL_To_Be_Added_
            for i in range(self.NUM_GMs-1):
                # noinspection PyUnresolvedReferences
                self._LOCAL_To_Be_Added_.append(self._LOCAL_To_Be_Added_[-1] + self.GV_LENs[i])

        self._freeze_self_()


    @property
    def NUM_GMs(self):
        """Return an int which represents how many GM this CGM has."""
        return self._NUM_GMs_

    @property
    def GLOBAL_num_dofs(self):
        """How many dofs in total in all cores."""
        return self._GLOBAL_num_dofs_

    @property
    def shape(self):
        """Return (local number of elements, local number of dofs)."""
        return len(self), sum(self.local_dofs_distribution)

    def __getitem__(self, item):
        """Return the chained gathering vector for the mesh element numbered by `item`.

        :param item: #item mesh element.
        :return: A 1d numpy.array. But we do not save it, it is made whenever we call this method.
        """
        if self.___NUM___ == 1:
            return self.GMs[0][item].full_vector
        else:
            _ = [gm[item].full_vector+self._To_Be_Added_[i] for i, gm in enumerate(self._GMs_)]
            return np.concatenate(_)

    def __iter__(self):
        """Go through all local mesh elements."""
        for i in self._GMs_[0]:
            yield i

    def __len__(self):
        """local length, equal to the amount of local elements."""
        return len(self._GMs_[0])

    def __contains__(self, item):
        """If a mesh element is contained by this CGM? Or if a mesh element is local?"""
        return item in self._GMs_[0]

    def __eq__(self, other):
        """If two CGMs are equal?"""
        if self is other:
            RETURN = True
        else:
            if other.__class__.__name__ != 'Chain_Gathering_Matrix':
                RETURN = False
            elif other.___NUM___ != self.___NUM___:
                RETURN = False
            else:
                _ = list()
                for i, gm in enumerate(self._GMs_):
                    tf = gm == other._GMs_[i]
                    _.append(tf)
                    if not tf:
                        break
                RETURN = all(_)

        RETURN = cOmm.allreduce(RETURN, op=MPI.LAND)
        return RETURN

    @property
    def mesh_type(self):
        return self._mesh_type_

    @property
    def GMs(self):
        return self._GMs_

    @property
    def local_ranges(self):
        if self._local_ranges_ is None:

            if len(self) == 0:
                self._local_ranges_ = [tuple() for _ in range(len(self.GMs))]

            else:
                self._local_ranges_ = list()
                tba = self._To_Be_Added_
                for i, gm in enumerate(self.GMs):
                    LR = gm.local_range
                    MIN, MAX = LR
                    MIN += tba[i]
                    MAX += tba[i]
                    self._local_ranges_.append((MIN, MAX))

        return self._local_ranges_


    @property
    def GLOBAL_dofs_distribution(self):
        """Like [100, 127, 80], we know GMs[0].GLOBAL_num_dofs=100, GMs[1].GLOBAL_num_dofs=127,
         GMs[2].GLOBAL_num_dofs=80. So we know, for example, how to distribution the result vector
         by solving the assembled EWC_SparseMatrix."""
        return self._dofs_distribution_


    @property
    def local_dofs_distribution(self):
        """Let local_dofs_distribution = [20, 12, 8], we then know
        GMs[0] has 20 local dofs, of shape (., 20), and
        GMs[1] has 12 local dofs, of shape (., 12) and so on.
        """
        return self._local_dofs_distribution_


    @property
    def local_dofs_ranges(self):
        """if self.local_dofs_distribution = [20, 12, 8], then self.local_dofs_ranges = [range(0,20), range(20, 32), range(32,40)]"""
        SUM = 0
        ranges = list()
        for i in self._local_dofs_distribution_:
            ranges.append(range(SUM, SUM+i))
            SUM += i
        return ranges


    @property
    def find(self):
        return self._FIND_

    @property
    def GV_LENs(self):
        """The length (size) of the gathering_vector (gv.full_vector) of each GM.

        Therefore, sum(GV_LENs) == len(self[i]).

        """
        if len(self) == 0:
            LENs = [None for _ in range(self.NUM_GMs)]
        else:
            for i in self:
                LENs = [len(gm[i]) for gm in self.GMs]
                break
        return LENs












if __name__ == '__main__':
    # mpiexec -n 6 python tools/linear_algebra/gathering/chain_matrix.py
    GV = Gathering_Vector
    GM = Gathering_Matrix