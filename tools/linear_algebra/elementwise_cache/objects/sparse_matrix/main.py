# -*- coding: utf-8 -*-
import types
from screws.freeze.main import FrozenOnly
from scipy import sparse as spspa
# from tools.linear_algebra.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix
# from tools.linear_algebra.gathering.irregular.ir_chain_matrix.main import iR_Chain_Gathering_Matrix
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.customize import SpaMat_Customize
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.adjust import SpaMat_Adjust
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.blocks.main import EWC_SpaMat_Blocks
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.condition.main import EWC_SpaMat_Condition

from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.matmul import ___MATMUL___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.vecmul import ___VECMUL___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.add import ___ADD___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.sub import ___SUB___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.truediv import ___TRUE_DIV___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.transpose import ___TRANSPOSE___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.inv import ___LinearAlgebraINV___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.neg import ___NEG___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.mul import ___MUL___
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.helpers.Psum import SpaMat_PRIVATE_sum

from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.assembler import EWC_SparseMatrix_Assembler
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.do import EWC_SparseMatrix_Do
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.IS import EWC_SparseMatrix_IS
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.visualize import EWC_SparseMatrix_Vis

from tools.linear_algebra.elementwise_cache.objects.column_vector.main import EWC_ColumnVector

from tools.linear_algebra.gathering.chain import GatheringMatrixChaining

class EWC_SparseMatrix(FrozenOnly):
    """
    Element-wise cached sparse matrix (2D).

    :param mesh_elements: If it is given as a mesh, we will get the elements from the mesh.
    :param data_generator:
        1) `data_generator = (int, int )` and `cache_key_generator = None `
            we make locally empty sparse matrix of shape `data_generator`.
        2)  `data_generator = (int, int )` and `cache_key_generator = constant `
            we make locally empty sparse matrix of shape `data_generator`. (just like situation 1).
        3) `data_generator = ('identity', int-a)` and `cache_key_generator = constant `
            We will make identity local sparse matrix of shape (int-a, int-a) in all mesh elements.
        4) `data_generator` is a dict. This means the data are already stored in this dict. We do
            not need to set `cache_key_generator` anymore of course.

    :param cache_key_generator:
        1) `cache_key_generator = 'all_diff'`
            The local sparse matrix will be all different in all mesh elements.
        2) `cache_key_generator = 'constant'`
            The local sparse matrix will be all same in all mesh elements.
        3) `cache_key_generator = 'no_cache'`
            We will not cache the sparse matrix.
        else: `cache_key_generator = else`:
            we have a customized `cache_key_generator`

        When `data_generator = (x, y)` where `x`, `y` are positive integers, we make it empty sparse
            matrix in all elements.
    :param bmat_shape: If this EWC instance is made from a bmat, `bmat_shape` will no longer be
        False, and it will become the bmat shape.

        Or, for example,

        if bmat_shape = False: it is not from a bmat
        else bmat_shape should be of shape (2,) and is representing the block shape. For example,
            M = bmat([[A, B, C], [D, E, None]]), then bmat_shape = [2,3].

    """
    def __init__(self, mesh_elements, data_generator, cache_key_generator=None, bmat_shape=False):
        """

        :param mesh_elements:
        :param data_generator:
        :type data_generator: list, tuple, callable
        :param cache_key_generator:
        :param bmat_shape: Do not set this. It is an indicator used to indicate if we are generating
            a EWC_matrix by `bmat` other EWC_matrices.
        """
        # check mesh elements ---------------------------------------------------------------
        if mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh_Elements', '_2dCSCG_Mesh_Elements'):
            self._elements_ = mesh_elements
        elif mesh_elements.__class__.__name__ in ('miUsGrid_TriangularMesh_Elements',):
            self._elements_ = mesh_elements
        elif mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'):
            self._elements_ = mesh_elements.elements
        elif mesh_elements.__class__.__name__ in ('miUsGrid_TriangularMesh',):
            self._elements_ = mesh_elements.elements
        elif isinstance(mesh_elements, dict):
            self._elements_ = mesh_elements
        else:
            raise Exception(f"{mesh_elements}")

        # we can accept a dictionary as a data generator, we will wrap it with a method -----------
        if isinstance(data_generator, dict):
            self.___fully_pre_data_DICT___ = True # the data are already created!

            assert len(data_generator) == len(self._elements_), "dict key wrong."
            for _ in data_generator: assert _ in self._elements_, "dict key wrong."
            self.___dict_DG___ = data_generator
            data_generator = self.___PRIVATE_dict_2_method_data_generator___

            if cache_key_generator is None:
                cache_key_generator = 'no_cache'
                # the data are in the dict anyway, so do not need to be cached.
        else:
            self.___fully_pre_data_DICT___ = False # the data are not created yet.

        #---------------parse data type ------------------------------------------------------------
        DATA_TYPE = None
        if isinstance(data_generator, (list, tuple)) and data_generator[0] == 'identity':
            # data_generator[1] = a (int), (a, a) be the shape of the local identity matrix.
            DATA_TYPE = "IDENTITY"

        elif isinstance(data_generator, (list, tuple)) and len(data_generator) == 2 and \
                all([data_generator[i] % 1 == 0 and data_generator[i] > 0 for i in range(2)]):
            # the `data_generator` is the shape of the empty local sparse matrix.
            DATA_TYPE = "EMPTY"
        else:
            pass

        #---- parse default cache_key_generator ----------------------------------------------------
        if cache_key_generator is None:

            if DATA_TYPE == 'IDENTITY':
                pass
            elif DATA_TYPE == 'EMPTY':
                pass
            else:
                cache_key_generator = 'constant'
        else:
            pass

        # ----  we are making identity sparse matrices ----------------------------------------
        if DATA_TYPE == "IDENTITY":

            SHAPE = data_generator[1]
            assert len(data_generator) == 2 and (SHAPE % 1 == 0 and SHAPE > 0), \
                f"`data_generator` = {data_generator} is wrong. To generate identity local matrix, " \
                f"use, for example, data_generator = ('identity, i) where i is an positive integer " \
                f"representing the shape, (i, i), of the local identity matrix." \

            self.___IDENTITY_SHAPE___ = SHAPE
            self._DG_ = self.___PRIVATE_identity_cache_data_generator___
            self._KG_ = self.___PRIVATE_constant_cache_key_generator___

        # we are making empty sparse matrices ---------------------------------------------------
        elif DATA_TYPE == "EMPTY":
            assert isinstance(data_generator, (list, tuple)) and len(data_generator) == 2, \
                f"When `cache_key_generator` is None, we make empty sparse matrix in all elements, thus " \
                f"`data_generator` must be a tuple or list of length 2. Now it is {data_generator}."

            assert all([data_generator[i] % 1 == 0 and data_generator[i] > 0 for i in range(2)]), \
                f"`data_generator` = {data_generator} is wrong. Two members should be int and > 0."

            self.___EMPTY_SHAPE___ = data_generator
            self._DG_ = self.___PRIVATE_empty_cache_data_generator___
            self._KG_ = self.___PRIVATE_constant_cache_key_generator___


        # we are making regular sparse matrices --------------------------------------------------
        elif DATA_TYPE is None: # regular
            if cache_key_generator == 'all_diff': # all elements return different things but still cache all.
                # although all different, we cache everything because it may be used over iterations.
                self._DG_ = data_generator
                self._KG_ = self.___PRIVATE_all_different_cache_key_generator___
            elif cache_key_generator == 'constant': # return the same sparse matrix for all elements.
                # the data_generator should be the data itself
                if spspa.isspmatrix_csc(data_generator) or spspa.isspmatrix_csr(data_generator):
                    self.___DGD___ = data_generator # save it, then we can call it.
                    self._DG_ = self.___PRIVATE_constant_cache_data_generator___
                else:
                    self._DG_ = data_generator
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

        else:
            raise NotImplementedError(f"cannot deal with data type = {DATA_TYPE}.")

        #--------------------------------------------------------------------------------------
        self._gathering_matrices_0_ = None
        self._gathering_matrices_1_ = None
        self._cache_ = dict() # do not use self.___PRIVATE_reset_cache___()
        self.___CT___ = '>CT<'
        self.___NC___ = '>NC<'
        self.___IS_CT___ = False
        self.___IS_NC___ = False
        self.____CT_DG____ = None # the cache for constant data.
        self.___CHECK_repeat_CT___ = True
        self.___CHECK_repeat_CT___ = True
        self.___repeat_CK___ = ''
        self._customize_ = SpaMat_Customize(self)
        self._bmat_shape_ = bmat_shape
        self._assembler_ = None
        self._do_ = None
        self._IS_ = None
        self._visualize_ = None
        self._adjust_ = None
        self._blocks_ = None
        self._condition_ = None
        self._freeze_self_()

    def __repr__(self):
        return f'EWC_SpaMat:{id(self)}'

    def ___PRIVATE_reset_cache___(self):
        self._cache_ = dict()
        self.assembler.___PRIVATE_reset_cache___()

    def ___PRIVATE_all_different_cache_key_generator___(self, i):
        """cache key will be different for all elements since we use their id."""
        return str(id(self._elements_[i]))

    # noinspection PyUnusedLocal
    def ___PRIVATE_constant_cache_key_generator___(self, i):
        return self.___CT___

    # noinspection PyUnusedLocal
    def ___PRIVATE_constant_cache_data_generator___(self, i):
        return self.___DGD___

    # noinspection PyUnusedLocal
    def ___PRIVATE_empty_cache_data_generator___(self, i):
        """"""
        return spspa.csr_matrix(self.___EMPTY_SHAPE___)

    # noinspection PyUnusedLocal
    def ___PRIVATE_identity_cache_data_generator___(self, i):
        return spspa.identity(self.___IDENTITY_SHAPE___, format='csr')

    # noinspection PyUnusedLocal
    def ___PRIVATE_no_cache_key_generator___(self, i):
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

        :param gathering_matrices:
            it can be like:
            1. gathering_matrices = (CGM1, CGM2) # CGM1, CGM2 be Chain_Gathering_Matrix.
            2. gathering_matrices = (u2, P3) # we will make two Chain_Gathering_Matrix from them.
            3. gathering_matrices = ([u2, P3], [u2, P3]) # we will make two Chain_Gathering_Matrix from them.
        :return:
        """

        GMS = list()
        for gm in gathering_matrices:
            if not isinstance(gm, (tuple, list)):
                name = gm.__class__.__name__

                if name in ('iR_Chain_Gathering_Matrix', 'Chain_Gathering_Matrix'):
                    GMS.append(gm)
                elif name in ('iR_Gathering_Matrix', 'Gathering_Matrix'):
                    GMS.append([gm,])
                elif hasattr(gm, '___IS_ADF___') and gm.___IS_ADF___:
                    GMS.append([gm.prime.numbering.gathering,])
                else:
                    GMS.append([gm.numbering.gathering,])
            else:
                GM = list()

                for _i in gm:
                    name = _i.__class__.__name__
                    if name in ('iR_Gathering_Matrix', 'Gathering_Matrix'):
                        GM.append(_i)
                    elif hasattr(_i, '___IS_ADF___') and _i.___IS_ADF___:
                        GM.append(_i.prime.numbering.gathering)
                    else:
                        GM.append(_i.numbering.gathering)

                GMS.append(GM)

        FINAL_GMS = list()
        for gms in GMS:

            if gms.__class__.__name__  == 'iR_Chain_Gathering_Matrix':

                assert gms.chain_method == self.assembler.chain_method, \
                    f"The GM.chain_method = {gms.chain_method} does not match that of the " \
                    f"vector, {self.assembler.chain_method}."

                FINAL_GMS.append(gms)

            elif gms.__class__.__name__ == 'Chain_Gathering_Matrix':

                assert gms.chain_method == self.assembler.chain_method, \
                    f"The GM.chain_method = {gms.chain_method} does not match that of the " \
                    f"vector, {self.assembler.chain_method}."

                FINAL_GMS.append(gms)

            elif isinstance(gms, list):

                FINAL_GMS.append(
                    GatheringMatrixChaining(*gms)(chain_method=self.assembler.chain_method)
                )

                # if all([_.__class__.__name__ == 'iR_Gathering_Matrix' for _ in gms]):
                #     FINAL_GMS.append(
                #         iR_Chain_Gathering_Matrix(gms, chain_method=self.assembler.chain_method)
                #     )
                #
                # elif all([_.__class__.__name__ == 'Gathering_Matrix' for _ in gms]):
                #     FINAL_GMS.append(
                #         Chain_Gathering_Matrix(gms, chain_method=self.assembler.chain_method)
                #     )
                #
                # else:
                #     raise Exception()
            else:
                raise Exception()

        self._gathering_matrices_0_, self._gathering_matrices_1_ = FINAL_GMS



    @property
    def IS(self):
        """The assembler"""
        if self._IS_ is None:
            self._IS_ = EWC_SparseMatrix_IS(self)
        return self._IS_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = EWC_SparseMatrix_Vis(self)
        return self._visualize_

    @property
    def adjust(self):
        if self._adjust_ is None:
            self._adjust_ = SpaMat_Adjust(self)
        return self._adjust_

    @property
    def blocks(self):
        if self._blocks_ is None:
            self._blocks_ = EWC_SpaMat_Blocks(self)
        return self._blocks_

    @property
    def do(self):
        """The assembler"""
        if self._do_ is None:
            self._do_ = EWC_SparseMatrix_Do(self)
        return self._do_

    @property
    def assembler(self):
        """The assembler"""
        if self._assembler_ is None:
            self._assembler_ = EWC_SparseMatrix_Assembler(self)
        return self._assembler_

    @property
    def condition(self):
        if self._condition_ is None:
            self._condition_ = EWC_SpaMat_Condition(self)
        return self._condition_

    @property
    def assembled(self):
        """
        We will call the assembler with the default routine to assemble self into a global matrix.

        :return:
        :rtype GlobalMatrix:
        """
        return self.assembler()

    def __len__(self):
        return len(self._elements_)

    def __contains__(self, item):
        return item in self._elements_

    def __iter__(self):
        for i in self._elements_:
            yield i

    def ___getitem_pre_customizing___(self, item):
        """"""
        assert item in self, "Out of range!"
        if self.___IS_NC___:
            RETURN = self._DG_(item)
        elif self.___IS_CT___:
            RETURN = self.____CT_DG____
        else:
            # noinspection PyCallingNonCallable
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
                # once reach here, we no longer do self._KG_(i) for further items because we know it is CT
            elif self.___NC___ in ck:
                # once it is or one component of it is not cached, we compute it every single time.
                RETURN = self._DG_(item)
                self.___IS_NC___ = True
                # once reach here, we no longer do self._KG_(i) for further items because we know it is NC
            else:
                if ck in self._cache_:
                    RETURN = self._cache_[ck]
                else:
                    RETURN = self._DG_(item)
                    self._cache_[ck] = RETURN

        return RETURN

    def __getitem__(self, item):
        """"""
        RETURN = self.___getitem_pre_customizing___(item)
        # customization is after the cache, so we can do whatever customization afterwards.
        RETURN = self.customize.___PRIVATE_do_execute_customization___(RETURN, item)
        return RETURN

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
        multiply other (int or float) with self, e.g. 7 * self.


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
        assert other.__class__.__name__ == 'EWC_SparseMatrix', f"other is a {other.__class__.__name__}"
        DKC = ___ADD___(self, other)
        return EWC_SparseMatrix(self._elements_, DKC.__DG_call__, DKC.__KG_call__)

    def __matmul__ (self, other):
        """"""
        if other.__class__.__name__ == 'EWC_SparseMatrix':
            DKC = ___MATMUL___(self, other)
            return EWC_SparseMatrix(self._elements_, DKC.__DG_call__, DKC.__KG_call__)


        elif other.__class__.__name__ == 'EWC_ColumnVector':
            DKC = ___VECMUL___(self, other)
            return EWC_ColumnVector(self._elements_, DKC.__DG_call__, DKC.__KG_call__)


        elif hasattr(other, 'standard_properties') and 'form' in other.standard_properties.tags:
            DKC = ___VECMUL___(self, other.cochain.EWC)
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
        """inv of self.

        Watch out, this could be very slow.
        """
        data_generator = ___LinearAlgebraINV___(self)
        return EWC_SparseMatrix(self._elements_, data_generator, self._KG_)

    def ___PRIVATE_sum___(self, others):
        """self + all others."""
        Mat = [self,] + others
        data_generator = SpaMat_PRIVATE_sum(Mat)
        # noinspection PyTypeChecker
        RETURN = EWC_SparseMatrix(self._elements_, data_generator, data_generator.__KG_call__)

        return RETURN