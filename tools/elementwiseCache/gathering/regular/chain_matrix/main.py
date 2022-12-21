# -*- coding: utf-8 -*-
""""""
import sys
if './' not in sys.path:
    sys.path.append('/')

from components.freeze.main import FrozenOnly
from root.config.main import *
from scipy.sparse import lil_matrix

from tools.elementwiseCache.gathering.vector import Gathering_Vector
from tools.elementwiseCache.gathering.regular.matrix.main import Gathering_Matrix
from tools.elementwiseCache.gathering.regular.chain_matrix.do.main import ___Chain_Gathering_Matrix_DO___


class Chain_Gathering_Matrix(FrozenOnly):
    """We chain some Gathering_Matrix together, for example, chain (GM0, GM1, GM2).

    Then in GM1, the numbering does not
    start with 0 anymore, it starts with ``GM0.GLOBAL_num_dofs``. As for GM2, its numbering now starts with
    ``GM0.GLOBAL_num_dofs + GM1.GLOBAL_num_dofs``.

    The idea is not just chain them, we also provide a mechanism to save memory, we do not make a new Gathering_Matrix
    for this chained Chain_Gathering_Matrix.
    """
    def __init__(self, GMs, chain_method=None):
        """
        GMs can be a gathering matrix or a chain gathering matrix. Or a list of them.

        :param GMs:
        """

        if GMs.__class__.__name__ in ('Gathering_Matrix', 'Chain_Gathering_Matrix,'):
            GMs = [GMs, ]

        assert isinstance(GMs, (list, tuple)) and len(GMs) >= 1, \
            "we need a list or tuple of at least one gathering matrices."

        mesh_types = list()
        NEW_GMs = list()
        for i, gm in enumerate(GMs):

            mesh_types.append(gm.mesh_type)
            assert mesh_types[i] == mesh_types[0], \
                f"mesh_types[{i}]: {mesh_types[i]} is different from mesh_types[0]: {mesh_types[0]}"

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

        To_Be_Added = [0, ]
        for gm in GMs[:-1]:
            GLOBAL_num_dofs = gm.global_num_dofs
            # noinspection PyUnresolvedReferences
            To_Be_Added.append(To_Be_Added[-1]+GLOBAL_num_dofs)

        self._To_Be_Added_ = To_Be_Added
        self._GLOBAL_num_dofs_ = To_Be_Added[-1] + GMs[-1].global_num_dofs

        self._GMs_ = GMs
        self._NUM_GMs_ = len(GMs)
        self._local_dofs_distribution_ = [_.global_shape[1] for _ in self.GMs]
        self._local_ranges_ = None
        self._do_ = None

        # -------- silly way is nice when len(GMs) == 1 ------------------------------------------------3
        if len(GMs) == 1 or chain_method is None:
            chain_method = 'silly'
        else:
            pass
        self._chain_method_ = chain_method

        # -----------------------------------------------------------------------------------------------3
        if chain_method == 'sequent':
            # ------------------- do the re-numbering ---------------------------------------------2
            for rank in range(SIZE):
                if RANK == rank:
                    # ---- initialize cache which saving the re-numbering ---------------------1
                    cache = list()
                    for gm in GMs:
                        cache.append(lil_matrix((1, gm.global_num_dofs), dtype=int))

                    # ------------- update cache ----------------------------------------------1
                    if RANK == 0:
                        current_number = 0
                        filters = [0 for _ in range(len(GMs))]
                    else:
                        current_number, filters = COMM.recv(source=RANK - 1, tag=RANK)

                    # ----- renumbering -------------------------------------------------------1
                    for i in self:
                        for j, gm in enumerate(GMs):
                            fv = gm[i].full_vector
                            Filter = filters[j]
                            fv = fv[fv >= Filter]  # these dofs should be re-numbered
                            amount_re_numbered = len(fv)
                            filters[j] += amount_re_numbered
                            cache[j][0, fv] = np.arange(current_number, current_number+amount_re_numbered)
                            current_number += amount_re_numbered

                    # --------------- send to next --------------------------------------------1
                    if RANK != SIZE - 1:
                        COMM.send([current_number, filters], dest=RANK + 1, tag=RANK + 1)
                    else:
                        pass

            # -------- merge cache ----------------------------------------------------------------2
            # noinspection PyUnboundLocalVariable
            cache = COMM.gather(cache, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                CACHE = list()
                for gm in GMs:
                    CACHE.append(lil_matrix((1, gm.global_num_dofs), dtype=int))

                for caches in cache:
                    for i, ch in enumerate(caches):
                        CACHE[i] += ch

                del caches, cache

            else:
                pass

            # ------------ distribute re-numbering to each core -----------------------------------2
            local_dofs = list()
            for gm in GMs:
                local_dofs.append(gm.local_dofs)

            local_dofs = COMM.gather(local_dofs, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                cache_d = list()
                for core in range(SIZE):
                    cache_d.append(
                        [lil_matrix((1, gm.global_num_dofs), dtype=int) for gm in GMs]
                    )

                for i, LDs in enumerate(local_dofs):
                    for j, lds in enumerate(LDs):
                        lds = list(lds)
                        # noinspection PyUnboundLocalVariable
                        cache_d[i][j][0, lds] = CACHE[j][0, lds]

            else:
                cache_d = None

            cache_d = COMM.scatter(cache_d, root=MASTER_RANK)

            del local_dofs
            self._cache_ = [_.tocsr() for _ in cache_d]

        # ------------------------------------------------------------------------------------------------3
        elif chain_method == 'silly':
            if len(self) > 0:
                for i in self:
                    assert sum(self.GV_LENs) == len(self[i]), "A safety check."
                    break
                self._LOCAL_To_Be_Added_ = [0, ]
                # when len(self) == 0, we do not even initialize self._LOCAL_To_Be_Added_
                for i in range(self.num_GMs - 1):
                    # noinspection PyUnresolvedReferences
                    self._LOCAL_To_Be_Added_.append(self._LOCAL_To_Be_Added_[-1] + self.GV_LENs[i])
            else:
                pass

        # -----------------------------------------------------------------------------------------------3
        else:
            raise NotImplementedError(f"chain_method={chain_method} not implemented")
        # ==============================================================================================3

        self._freeze_self_()

    @property
    def chain_method(self):
        """Which method we use to chain the gathering matrices."""
        return self._chain_method_

    @property
    def ___Pr_IS_regular___(self):
        return True

    @property
    def num_GMs(self):
        """Return an int which represents how many GM this CGM has."""
        return self._NUM_GMs_

    @property
    def global_num_dofs(self):
        """How many dofs in total in all cores."""
        return self._GLOBAL_num_dofs_

    def __getitem__(self, item):
        """Return the chained gathering vector for the mesh element numbered by `i`.

        :param item: #i mesh element.
        :return: A 1d numpy.array. But we do not save it, it is made whenever we call this method.
        """
        if self.___NUM___ == 1:
            return self.GMs[0][item].full_vector

        else:
            if self.chain_method == 'silly':
                _ = [gm[item].full_vector+self._To_Be_Added_[i] for i, gm in enumerate(self._GMs_)]
                return np.concatenate(_)

            elif self.chain_method == 'sequent':
                _ = list()
                for j, gm in enumerate(self.GMs):
                    old_numbering = gm[item].full_vector
                    new_numbering = self._cache_[j][0, old_numbering].toarray().ravel()
                    _.append(new_numbering)
                return np.concatenate(_)

            else:
                raise NotImplementedError(f"chain method={self.chain_method} is not implemented.")

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
                RETURN = True
                for i, gm in enumerate(self._GMs_):
                    tf = gm == other._GMs_[i]
                    if not tf:
                        RETURN = False
                        break

        RETURN = COMM.allreduce(RETURN, op=MPI.LAND)
        return RETURN

    @property
    def mesh_type(self):
        return self._mesh_type_

    @property
    def GMs(self):
        return self._GMs_

    @property
    def local_ranges(self):
        """The local dofs are in this range. This only make sense for silly chain method, try to
        avoid this in the code.
        """
        if self.chain_method == 'silly':
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
        else:
            raise Exception(f'local_ranges only applicable to silly chain_matrix.')


    @property
    def ___Pr_regular__local_dofs_ranges___(self):
        """if self._local_dofs_distribution_ = [20, 12, 8],
        then self.___Pr_regular__local_dofs_ranges___ = [range(0,20), range(20, 32), range(32,40)]"""
        SUM = 0
        ranges = list()
        for i in self._local_dofs_distribution_:
            ranges.append(range(SUM, SUM+i))
            SUM += i
        return ranges


    @property
    def do(self):
        if self._do_ is None:
            self._do_ = ___Chain_Gathering_Matrix_DO___(self)
        return self._do_

    @property
    def GV_LENs(self):
        """The length (size) of the gathering_vector (gv.full_vector) of each GM.

        So if GV_LENs = [10, 15, 8], we know that the frst 10 numbers are the numbering of the
        GMs[0], GMs[1] has 15 dofs, and  GMs[2] has 8 dofs.

        Therefore, sum(GV_LENs) == len(self[i]).

        """
        LENs = [None for _ in range(self.num_GMs)]
        if len(self) == 0:
            pass
        else:
            for i in self:
                LENs = [len(gm[i]) for gm in self.GMs]
                break
        return LENs


if __name__ == '__main__':
    # mpiexec -n 6 python tools/linear_algebra/gathering/chain_matrix.py
    GV = Gathering_Vector
    GM = Gathering_Matrix
