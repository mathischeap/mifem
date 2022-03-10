
import types
from screws.freeze.main import FrozenOnly
from scipy import sparse as spspa
from scipy.sparse import linalg as spspalinalg
from root.config.main import *
from tools.linear_algebra.data_structures.global_matrix.main import GlobalMatrix
from tools.linear_algebra.gathering.chain_matrix.main import Chain_Gathering_Matrix
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.customize import SpaMat_Customize

from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector



class EWC_SparseMatrix(FrozenOnly):
    """
    Element-wise cached sparse matrix (2D).

    :param mesh_elements:
    :param data_generator: When `data_generator = (x, y)` where `x`, `y` are positive integers, we make it empty sparse
        matrix in all elements.
    :param cache_key_generator:
    :param bmat_shape: If this EWC instance is made from a bmat and if yes, what is the bmat shape.
        if bmat_shape = False: it is not from a bmat
        else bmat_shape should be of shape (2,) and is representing the block shape. For example,
            M = bmat([[A, B, C], [D, E, None]]), then bmat_shape = [2,3].
    :type bmat_shape: list, tuple

    """
    def __init__(self, mesh_elements, data_generator, cache_key_generator=None, bmat_shape=False):

        # check mesh elements ---------------------------------------------------------------
        if mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh_Elements', '_2dCSCG_Mesh_Elements'):
            self._elements_ = mesh_elements
        elif mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'):
            self._elements_ = mesh_elements.elements
        else:
            raise Exception()
        self.____CT_DG____ = None

        # we can accept a dictionary as a data generator, we will wrap it with a method ...
        if isinstance(data_generator, dict):
            self.___dict_DG___ = data_generator
            data_generator = self.___PRIVATE_dict_2_method_data_generator___
        else:
            pass

        # Check if we make empty sparse matrices ------------------------------------------------
        if cache_key_generator is None:# We make an empty sparse matrix of shape `data_generator` in each element.
            EMPTY = True
        elif cache_key_generator == 'constant' and isinstance(data_generator, (list, tuple)):
            if len(data_generator) == 2 and \
                all([data_generator[i] % 1 == 0 and data_generator[i] > 0 for i in range(2)]):
                EMPTY = True
            else:
                EMPTY = False
        else:
            EMPTY = False

        # If we are making empty sparse matrices ---------------------------------------------------
        if EMPTY:
            assert isinstance(data_generator, (list, tuple)) and len(data_generator) == 2, \
                f"When `cache_key_generator` is None, we make empty sparse matrix in all elements, thus " \
                f"`data_generator` must be a tuple or list of length 2. Now it is {data_generator}."

            assert all([data_generator[i] % 1 == 0 and data_generator[i]>0 for i in range(2)]), \
                f"`data_generator` = {data_generator} is wrong. Two members should be int and > 0."

            self.___EMPTY_SHAPE___ = data_generator
            self._DG_ = self.___PRIVATE_empty_cache_data_generator___
            self._KG_ = self.___PRIVATE_constant_cache_key_generator___

        # we are not making empty sparse matrices --------------------------------------------------
        else:
            if cache_key_generator == 'all_diff': # all elements return different things but still cache all.
                # although all different, we cache everything because it may be used over iterations.
                self._DG_ = data_generator
                self._KG_ = self.___PRIVATE_all_different_cache_key_generator___
            elif cache_key_generator == 'constant': # return the same sparse matrix for all elements.
                # the data_generator should be the data itself
                assert spspa.isspmatrix_csc(data_generator) or spspa.isspmatrix_csr(data_generator)
                # must be a sparse csc or csr matrix.
                self.___DGD___ = data_generator # save it, then we can call it.
                self._DG_ = self.___PRIVATE_constant_cache_data_generator___
                self._KG_ = self.___PRIVATE_constant_cache_key_generator___
            elif cache_key_generator == 'no_cache': # do not cache for any elements.
                # use this when nothing is the same in elements and iterations: i.e. for the cross product
                self._DG_ = data_generator
                self._KG_ = self.___PRIVATE_no_cache_key_generator___
            else:
                # if reach here, cache_key_generator must be a method have one input (#element) apart from self.
                if isinstance(cache_key_generator, types.MethodType):
                    # noinspection PyUnresolvedReferences
                    assert cache_key_generator.__code__.co_argcount == 2
                else:
                    pass
                self._DG_ = data_generator
                self._KG_ = cache_key_generator

        self._gathering_matrices_0_ = None
        self._gathering_matrices_1_ = None
        self.___PRIVATE_reset_cache___()
        self.___NC___ = '>NC<'
        self.___CT___ = '>CT<'
        self.___IS_NC___ = False
        self.___IS_CT___ = False
        self.___CHECK_repeat_CT___ = True
        self.___repeat_CK___ = ''
        self._customize_ = SpaMat_Customize(self)
        self._bmat_shape_ = bmat_shape
        self._shape_ = None
        self._freeze_self_()


    def ___PRIVATE_reset_cache___(self):
        self._cache_ = dict()



    def ___PRIVATE_all_different_cache_key_generator___(self, i):
        """cache key will be different for all elements since we use their id."""
        return str(id(self._elements_[i]))

    def ___PRIVATE_constant_cache_key_generator___(self, i):
        assert i in self.elements
        return self.___CT___

    def ___PRIVATE_constant_cache_data_generator___(self, i):
        assert i in self.elements
        return self.___DGD___

    def ___PRIVATE_empty_cache_data_generator___(self, i):
        """"""
        assert i in self.elements
        return spspa.csr_matrix(self.___EMPTY_SHAPE___)

    def ___PRIVATE_no_cache_key_generator___(self, i):
        assert i in self.elements
        return self.___NC___

    def ___PRIVATE_dict_2_method_data_generator___(self, i):
        """When we get a dict as the data generator, we make wrap it with a method."""
        return self.___dict_DG___[i]


    @property
    def elements(self):
        """The mesh elements."""
        return self._elements_

    @property
    def gathering_matrices(self):
        """Return two Chain_Gathering_Matrix represent the row direction and column direction."""
        return self._gathering_matrices_0_, self._gathering_matrices_1_

    @gathering_matrices.setter
    def gathering_matrices(self, gathering_matrices):
        """Two Chain_Gathering_Matrix. If not, we make them to be Chain_Gathering_Matrix.

        :param gathering_matrices: it can be like:
            1. gathering_matrices = (CGM1, CGM2) # CGM1, CGM2 be Chain_Gathering_Matrix.
            2. gathering_matrices = (u2, P3) # we will make two Chain_Gathering_Matrix from them.
            3. gathering_matrices = ([u2, P3], [u2, P3]) # we will make two Chain_Gathering_Matrix from them.
        :return:
        """
        CGM0, CGM1 = gathering_matrices
        if CGM0.__class__.__name__ == 'Chain_Gathering_Matrix':
            pass
        else:
            if not isinstance(CGM0, (list, tuple)):
                CGM0 = [CGM0,]
            cgm0 = list()
            for _ in CGM0:
                if _.__class__.__name__ == 'Gathering_Matrix':
                    cgm0.append(_)
                else:
                    cgm0.append(_.numbering.gathering)
            CGM0 = Chain_Gathering_Matrix(cgm0)

        if CGM1.__class__.__name__ == 'Chain_Gathering_Matrix':
            pass
        else:
            if not isinstance(CGM1, (list, tuple)):
                CGM1 = [CGM1,]
            cgm1 = list()
            for _ in CGM1:
                if _.__class__.__name__ == 'Gathering_Matrix':
                    cgm1.append(_)
                else:
                    cgm1.append(_.numbering.gathering)
            CGM1 = Chain_Gathering_Matrix(cgm1)

        assert CGM0.__class__.__name__ == 'Chain_Gathering_Matrix', "I need Chain_Gathering_Matrix!"
        assert CGM1.__class__.__name__ == 'Chain_Gathering_Matrix', "I need Chain_Gathering_Matrix!"

        self._gathering_matrices_0_ = CGM0
        self._gathering_matrices_1_ = CGM1

        if self._shape_ is not None:
            assert self._gathering_matrices_0_.shape + \
                   self._gathering_matrices_1_.shape[1:] == self._shape_


    @property
    def assembled(self):
        """
        Assemble self in a global matrix.

        :return:
        :rtype GlobalMatrix:
        """
        assert self.gathering_matrices[0] is not None, "I have no gathering matrix"
        assert self.gathering_matrices[1] is not None, "I have no gathering matrix"
        GI = self.gathering_matrices[0]
        GJ = self.gathering_matrices[1]
        DEP = GI.GLOBAL_num_dofs
        WID = GJ.GLOBAL_num_dofs
        ROW = list()
        COL = list()
        DAT = list()

        A = spspa.csc_matrix((DEP, WID)) # initialize a sparse matrix

        for i in self:
            Mi = self[i]
            indices = Mi.indices
            indptr = Mi.indptr
            data = Mi.data
            nums = np.diff(indptr)
            if Mi.__class__.__name__ == 'csc_matrix':
                for j, num in enumerate(nums):
                    idx = indices[indptr[j]:indptr[j+1]]
                    ROW.extend(GI[i][idx])
                    COL.extend([GJ[i][j],]*num)
            elif Mi.__class__.__name__ == 'csr_matrix':
                for j, num in enumerate(nums):
                    idx = indices[indptr[j]:indptr[j+1]]
                    ROW.extend([GI[i][j],]*num)
                    COL.extend(GJ[i][idx])
            else:
                raise Exception("I can not handle %r."%Mi)
            DAT.extend(data)

            if len(DAT) > 1e7: # every 5 million data, we make it into sparse matrix.
                _ = spspa.csc_matrix((DAT, (ROW, COL)), shape=(DEP, WID)) # we make it into sparse
                del ROW, COL, DAT
                A += _
                del _
                ROW = list()
                COL = list()
                DAT = list()

        _ = spspa.csc_matrix((DAT, (ROW, COL)), shape=(DEP, WID))  # we make it into sparse
        del ROW, COL, DAT
        A += _
        del _

        assert A.shape == (GI.GLOBAL_num_dofs, GJ.GLOBAL_num_dofs)

        return GlobalMatrix(A)

    @property
    def GLOBAL_len(self):
        return self.elements.GLOBAL_num

    def __len__(self):
        return len(self._elements_)

    def __contains__(self, item):
        return item in self._elements_

    def __iter__(self):
        for i in self._elements_:
            yield i

    def __getitem__(self, item):
        """"""
        assert item in self, "Out of range!"
        if self.___IS_NC___:
            RETURN = self._DG_(item)
        elif self.___IS_CT___:
            RETURN = self.____CT_DG____
        else:
            ck = self._KG_(item)

            if self.___CHECK_repeat_CT___:
                if self.___CT___ in ck:
                    # if ck = '>CT<>CT<...', repeat_CK will be '>CT<'
                    temp = (ck + ck).find(ck, 1, -1)
                    if temp != -1:
                        self.___repeat_CK___ = ck[:temp]
                    # ...
                self.___CHECK_repeat_CT___ = False # only do above check once.

            if ck == self.___CT___ or self.___repeat_CK___ == self.___CT___:
                assert self.____CT_DG____ is None, "self.____CT_DG____ must be None so far"
                # one more cache to make it always cached even after operators
                self.____CT_DG____ = self._DG_(item)
                RETURN = self.____CT_DG____ # then we do not call the data generator
                self.___IS_CT___ = True
                # once reach here, we no longer do self._KG_(item) for further items because we know it is CT
            elif self.___NC___ in ck:
                # once it is or one component of it is not cached, we compute it every single time.
                RETURN = self._DG_(item)
                self.___IS_NC___ = True
                # once reach here, we no longer do self._KG_(item) for further items because we know it is NC
            else:
                if ck in self._cache_:
                    RETURN = self._cache_[ck]
                else:
                    RETURN = self._DG_(item)
                    self._cache_[ck] = RETURN

        # customization is after the cache, so we can do whatever customization afterwards.
        RETURN = self.customize.___PRIVATE_do_execute_customization___(RETURN, item)

        # here we do not check the shape of RETURN because we could use this method to compute shape, so if we do, we may end up in dead loop.

        return RETURN

    @property
    def shape(self):
        """The local shape : == (len(self),) + np.shape(self[i]). So the
        first value refer to how many local mesh elements. The second
        value refers to how many rows in the local sparse matrix. The third
        one refers to how many cols in the local sparse matrix.

        :return: A tuple of 3 integers.
        """
        if self._shape_ is not None: return self._shape_

        if self._gathering_matrices_0_ is not None and \
            self._gathering_matrices_1_ is not None:
            self._shape_ = self._gathering_matrices_0_.shape + \
                           self._gathering_matrices_1_.shape[1:]
        else:
            shape = None
            for i in self:
                Vi = self[i]
                shape = np.shape(Vi)
                break

            shape = cOmm.gather(shape, root=mAster_rank)
            _s_ = None
            if rAnk == mAster_rank:
                for s in shape:
                    if s is None:
                        pass
                    else:
                        if _s_ is None:
                            _s_ = s
                        else:
                            assert _s_ == s
                assert _s_ is not None, f"we must have found a _s_."
            else:
                pass

            _s_ = cOmm.bcast(_s_, root=mAster_rank)

            self._shape_ = (len(self),) + _s_

        return self._shape_


    @property
    def customizations(self):
        """All the customizations that have been added to me."""
        return self.customize._customizations_

    @property
    def customize(self):
        """We use sub-methods of these properties to add customization. These customizations will be
        executed when I am called."""
        return self._customize_

    @property
    def bmat_shape(self):
        return self._bmat_shape_

    def __mul__(self, other):
        """
        multiply self with other int of float, a * 7.

        :param other:
        :return:
        """
        data_generator = ___MUL___(self, other)
        RETURN = EWC_SparseMatrix(self._elements_, data_generator, self._KG_)
        if self.gathering_matrices != (None, None):
            RETURN.gathering_matrices = self.gathering_matrices
        return RETURN

    def __rmul__(self, other):
        """
        multiply other (int or float) with self, e.g. 7 * a.


        :param other:
        :return:
        """
        data_generator = ___MUL___(self, other)
        RETURN = EWC_SparseMatrix(self._elements_, data_generator, self._KG_)
        if self.gathering_matrices != (None, None):
            RETURN.gathering_matrices = self.gathering_matrices
        return RETURN

    def __truediv__(self, other):
        """
        division by int/float, like, a / 7.

        :param other:
        :return:
        """
        data_generator = ___TRUE_DIV___(self, other)
        RETURN = EWC_SparseMatrix(self._elements_, data_generator, self._KG_)
        if self.gathering_matrices != (None, None):
            RETURN.gathering_matrices = self.gathering_matrices
        return RETURN

    def __sub__(self, other):
        """self - EWC_SparseMatrix"""
        assert other.__class__.__name__ == 'EWC_SparseMatrix'
        assert self._elements_._mesh_ == other._elements_._mesh_
        DKC = ___SUB___(self, other)
        return EWC_SparseMatrix(self._elements_, DKC.__DG_call__, DKC.__KG_call__)

    def __neg__(self):
        """- EWC_SparseMatrix"""
        data_generator = ___NEG___(self)
        RETURN = EWC_SparseMatrix(self._elements_, data_generator, self._KG_)
        if self.gathering_matrices != (None, None):
            RETURN.gathering_matrices = self.gathering_matrices
        return RETURN

    def __add__(self, other):
        """self + EWC_SparseMatrix"""
        assert other.__class__.__name__ == 'EWC_SparseMatrix'
        assert self._elements_._mesh_ == other._elements_._mesh_
        DKC = ___ADD___(self, other)
        return EWC_SparseMatrix(self._elements_, DKC.__DG_call__, DKC.__KG_call__)

    def __matmul__ (self, other):
        if other.__class__.__name__ == 'EWC_SparseMatrix':
            assert self._elements_._mesh_ == other._elements_._mesh_
            DKC = ___MATMUL___(self, other)
            return EWC_SparseMatrix(self._elements_, DKC.__DG_call__, DKC.__KG_call__)
        elif other.__class__.__name__ == 'EWC_ColumnVector':
            DKC = ___VECMUL___(self, other)
            return EWC_ColumnVector(self._elements_, DKC.__DG_call__, DKC.__KG_call__)
        else:
            raise NotImplementedError()

    @property
    def T(self):
        """Transpose of self."""
        data_generator = ___TRANSPOSE___(self)
        RETURN = EWC_SparseMatrix(self._elements_, data_generator, self._KG_)

        if self.gathering_matrices != (None, None):
            RETURN.gathering_matrices = (self.gathering_matrices[1], self.gathering_matrices[0])

        return RETURN

    @property
    def inv(self):
        """inv of self."""
        data_generator = ___LinearAlgebraINV___(self)
        return EWC_SparseMatrix(self._elements_, data_generator, self._KG_)









class ___MUL___(FrozenOnly):
    def __init__(self, ewc, number):
        assert isinstance(number, (int, float))
        self._ewc_ = ewc
        self._number_ = number
        self._freeze_self_()

    def __call__(self, item):
        return self._ewc_[item] * self._number_

class ___TRUE_DIV___(FrozenOnly):
    def __init__(self, ewc, number):
        self._ewc_ = ewc
        self._number_ = number
        self._freeze_self_()

    def __call__(self, item):
        return self._ewc_[item] / self._number_

class ___SUB___(FrozenOnly):
    def __init__(self, EWC1, EWC2):
        self._ewc1_ = EWC1
        self._ewc2_ = EWC2
        self._freeze_self_()

    def __DG_call__(self, item):
        return self._ewc1_[item] - self._ewc2_[item]

    def __KG_call__(self, item):
        return self._ewc1_._KG_(item) + self._ewc2_._KG_(item)

class ___NEG___(FrozenOnly):
    def __init__(self, ewc):
        self._ewc_ = ewc
        self._freeze_self_()

    def __call__(self, item):
        return - self._ewc_[item]

class ___TRANSPOSE___(FrozenOnly):
    def __init__(self, ewc):
        self._ewc_ = ewc
        self._freeze_self_()

    def __call__(self, item):
        return self._ewc_[item].T

class ___LinearAlgebraINV___(FrozenOnly):
    def __init__(self, ewc):
        self._ewc_ = ewc
        self._freeze_self_()

    def __call__(self, item):
        return spspalinalg.inv(self._ewc_[item])

class ___ADD___(FrozenOnly):
    def __init__(self, EWC1, EWC2):
        self._ewc1_ = EWC1
        self._ewc2_ = EWC2
        self._freeze_self_()

    def __DG_call__(self, item):
        return self._ewc1_[item] + self._ewc2_[item]

    def __KG_call__(self, item):
        return self._ewc1_._KG_(item) + self._ewc2_._KG_(item)

class ___MATMUL___(FrozenOnly):
    def __init__(self, EWC1, EWC2):
        self._ewc1_ = EWC1
        self._ewc2_ = EWC2
        self._freeze_self_()

    def __DG_call__(self, item):
        return self._ewc1_[item] @ self._ewc2_[item]

    def __KG_call__(self, item):
        return self._ewc1_._KG_(item) + self._ewc2_._KG_(item)

class ___VECMUL___(FrozenOnly):
    def __init__(self, EWC_S, EWC_V):
        self._ewc_S_ = EWC_S
        self._ewc_V_ = EWC_V
        self._freeze_self_()

    def __DG_call__(self, item):
        return self._ewc_S_[item] @ self._ewc_V_[item]

    def __KG_call__(self, item):
        return self._ewc_S_._KG_(item) + self._ewc_V_._KG_(item)