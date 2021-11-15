# -*- coding: utf-8 -*-
""""""
from SCREWS.frozen import FrozenOnly
from TOOLS.__DEPRECATED__.assembler.main import GatheringMatrix
from itertools import chain
from root.config import *
from functools import lru_cache
from SCREWS.decorators import accepts


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
        """Return the Chain_Gathering_Vector of element #item.

        :param item: #item mesh element.
        :return: A 1d numpy.array. But we do not save it, it is made whenever we call this method.
        """
        if self.___NUM___ == 1:
            return self.GMs[0][item].full_vector
        else:
            _ = [gm[item].full_vector+self._To_Be_Added_[i] for i, gm in enumerate(self._GMs_)]
            return np.concatenate(_)

    def __iter__(self):
        for i in self._GMs_[0]:
            yield i

    def __len__(self):
        """local length, equal to the amount of local elements."""
        return len(self._GMs_[0])

    def __contains__(self, item):
        return item in self._GMs_[0]

    def __eq__(self, other):
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
        """The dofs in this core is in this range."""
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
    def FIND(self):
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


class ___Chain_Gathering_Matrix_FIND___(FrozenOnly):
    """"""
    def __init__(self, CGM):
        self._CGM_ = CGM
        self._freeze_self_()

    @accepts('self', (int, float, 'int32', 'int64'))
    def elements_contain_dof_numbered(self, m, N=None):
        """Find the element which contains the dof numbered `m`.

        :param m: the dof to be found.
        :param N: The max number of elements to find. If we already find `N` elements containing the dof `m`, we stop
            finding. When `N=None`, then we always check all local elements.
        :return: Return a list of ints or None

            An list of int(s) represent the local element number(s) that has the dof numbered `m` or None if no
            local element  contains that dof. Therefore, maybe in multiple cores, the method does not return None.

            We also return `where` (int): GMs[where] has the dof `m`. When the first output is None, `where = -1`.
        """
        assert (m % 1) == 0, f"m={m} is wrong."
        assert -self._CGM_.GLOBAL_num_dofs <= m < self._CGM_.GLOBAL_num_dofs, \
                f"dof numbered = {m} is out of range, it should be in " \
                f"[{-self._CGM_.GLOBAL_num_dofs}, {self._CGM_.GLOBAL_num_dofs-1}] " \
                f"cause the maximum global numbering for this chain_gathering_matrix " \
                f"is {self._CGM_.GLOBAL_num_dofs-1}."

        if m < 0: m += self._CGM_.GLOBAL_num_dofs
        assert 0 <= m < self._CGM_.GLOBAL_num_dofs, f"m={m} is wrong."
        if not isinstance(m, int): m = int(m)

        if N is None:
            if self._CGM_.mesh_type == '_2dCSCG': # if in a _3dCSCG mesh, of course, dofs can mostly shared by 4 elements.
                N = 4
            elif self._CGM_.mesh_type == '_3dCSCG': # if in a _3dCSCG mesh, of course, dofs can mostly shared by 8 elements.
                N = 8
            else:
                pass
        else:
            pass

        assert N is None or N > 0, f"At least we search for 1 element, right? Now it is N={N}"

        if len(self._CGM_) == 0: # chain_gathering_matrix is empty: no local element at all.
            return None

        LRS = self._CGM_.local_ranges

        where = -1

        if self._CGM_.NUM_GMs == 1:
            MIN, MAX = LRS[0]
            if MIN <= m < MAX:
                where = 0
            else:
                return None

        else:
            for i, LR in enumerate(LRS):
                MIN, MAX = LR
                if i == 0: # the first GM
                    if m < MIN:
                        return None
                    elif MIN <= m < MAX:
                        where = i
                        break
                    else:
                        pass

                else:

                    if LRS[i-1][1] <= m < MIN:
                        return None
                    elif MIN <= m < MAX:
                        where = i
                        break
                    else:
                        if i == self._CGM_.NUM_GMs - 1: # the last GM
                            if m >= MAX:
                                return None
                        else:
                            pass

        # if we reach here, we know this core may have element(s) containing the target dof...
        assert 0 <= where < self._CGM_.NUM_GMs, "MUST BE!" # we know we only need to look at GMs[where] to find m.

        m -= self._CGM_._To_Be_Added_[where]
        GM = self._CGM_.GMs[where]

        # now we go through `GM` to find the element(s) having `m`
        ELE = list()

        if N is None:
            for i in GM:
                if m in GM[i]: # NICE! we find one element containing m.
                    ELE.append(int(i))

        else:
            n = 0
            for i in GM:
                if m in GM[i]: # NICE! we find one element containing m.
                    ELE.append(int(i))
                    n += 1

                    if n == N: # we have find enough elements. Lets break the loop.
                        break

        if ELE == list():
            return None
        else:
            return ELE, int(where)


    @accepts('self', (int, float, 'int32', 'int64'))
    @lru_cache(maxsize=256)
    def elements_and_local_indices_of_dof(self, m, N=None):
        """

        :param m:
        :param N:
        :return: If no element contains dof m, return None. Else return two lists, one is for element(s), one is
            for the corresponding local index(s).
        """
        assert (m % 1) == 0, f"m={m} is wrong."
        assert -self._CGM_.GLOBAL_num_dofs <= m < self._CGM_.GLOBAL_num_dofs, \
                f"dof numbered = {m} is out of range, it should be in " \
                f"[{-self._CGM_.GLOBAL_num_dofs}, {self._CGM_.GLOBAL_num_dofs-1}] " \
                f"cause the maximum global numbering for this chain_gathering_matrix " \
                f"is {self._CGM_.GLOBAL_num_dofs-1}."
        if m < 0: m += self._CGM_.GLOBAL_num_dofs
        assert 0 <= m < self._CGM_.GLOBAL_num_dofs, f"m={m} is wrong."

        if not isinstance(m, int): m = int(m)

        OUT = self.elements_contain_dof_numbered(m, N=N)
        if OUT is None:
            return None
        elements, where = OUT
        local_indices = list()

        assert isinstance(elements, list), "A trivial check."

        m -= self._CGM_._To_Be_Added_[where]
        GM = self._CGM_.GMs[where]
        for e in elements:

            index = np.argwhere(GM[e].full_vector==m)[0,0] + self._CGM_._LOCAL_To_Be_Added_[where]
            local_indices.append(int(index))

        return elements, local_indices




















class Gathering_Matrix(GatheringMatrix):
    """A gathering matrix is a bunch of gathering vectors.

    Gathering vectors (as values) stored in a dictionary whose keys are mesh element indices.
    So it can not be used as, for example, trace-element-wise gathering matrix.

    :param dict gvd: A dictionary of gathering vectors.
    """
    def __init__(self, gvd, mesh_type=None):
        assert isinstance(gvd, dict), "I need a dictionary."
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
        self._freeze_self_()

    def __getitem__(self, item):
        return self._gvd_[item]

    def __iter__(self):
        for i in self._gvd_:
            yield i

    def __contains__(self, item):
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
        """Return (global number of elements, local number of dofs)."""
        return self._mesh_type_

    @property
    def local_range(self):
        """The range of numbering of dofs in this core."""
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

    def DO_hstack(self, *args):
        """
        Hstack other gathering matrices to the right of self.

        For example, we can do GM1.DO_hstack(GM2, GM3).

        :param args:
        :return:
        """
        for ogm in args:
            assert ogm.__class__.__name__ == 'Gathering_Matrix', f"I am stacking a {ogm.__class__.__name__}, wrong!"
            assert len(self) == len(ogm), "Length dis-match."
        for i in self:
            for ogm in args:
                assert i in ogm, "elements dis-match."

        ____ = [self.GLOBAL_num_dofs, *[ogm.GLOBAL_num_dofs for ogm in args]][:-1]
        upon = list()
        for j, _ in enumerate(args):
            upon.append(sum(____[0:j+1]))

        gvd = dict()
        for i in self:
            gvd_i_gv = self[i].full_vector
            for j, ogm in enumerate(args):
                gvd_i_gv = np.concatenate([gvd_i_gv, ogm[i].full_vector + upon[j]])
            gvd[i] = Gathering_Vector(i, gvd_i_gv)
        return Gathering_Matrix(gvd)









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
            self._gv_ = np.array([i for i in CHAIN])
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
        if self._full_vector_ is None:
            if self._gv_.__class__.__name__ == 'ndarray':
                self._full_vector_ = self._gv_
            elif isinstance(self._gv_, range):
                self._full_vector_ = np.array([j for j in self])
            else:
                raise Exception()
        return self._full_vector_

    def __getitem__(self, item):
        return self.full_vector[item]

    def __contains__(self, item):
        return item in self._gv_

    def __iter__(self):
        for j in self._gv_:
            yield j

    def __len__(self):
        return len(self._gv_)

    def __eq__(self, other):
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

    def ___PRIVATE_find_max_label___(self):
        return np.max(self.full_vector)

    def ___PRIVATE_find_min_label___(self):
        return np.min(self.full_vector)


