

import types
from SCREWS.frozen import FrozenOnly
from scipy import sparse as spspa
from scipy.sparse import linalg as spspalinalg
from root.config import *
from TOOLS.linear_algebra.data_structures import GlobalMatrix, GlobalVector
from TOOLS.linear_algebra.gathering import Chain_Gathering_Matrix
from SCREWS.decorators import accepts


class EWC_ColumnVector(FrozenOnly):
    """
    Element-wise cached column vector (csc_matrix of shape (n,1)).

    example::

        EWC_ColumnVector(mesh, 5) : generate an empty vector of shape (5,1) for all elements.


    :param mesh_elements:
    :param data_generator:
    :param con_shape:
    :type con_shape: bool, tuple, list
    """
    def __init__(self, mesh_elements, data_generator, cache_key_generator=None, con_shape=False):
        if mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh_Elements', '_2dCSCG_Mesh_Elements'):
            self._elements_ = mesh_elements
        elif mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'):
            self._elements_ = mesh_elements.elements
        else:
            raise Exception()
        self.____CT_DG____ = None

        EMPTY = False
        IS_Number = False

        if data_generator.__class__.__name__ in ('int', 'float', 'int32', 'int64'):
            IS_Number = True

        if cache_key_generator is None and IS_Number: # make empty sparse vector.
            # noinspection PyTypeChecker
            if data_generator % 1 == 0:

                EMPTY = True
                self.___EMPTY_LENGTH___ = data_generator

            else:
                raise Exception(f"empty sparse vector of local length={data_generator} is wrong.")

        elif cache_key_generator is None:
            cache_key_generator = 'no_cache' # DO NOT CHANGE THIS DEFAULT SETTING!!!

        else:
            pass

        if EMPTY:
            self._DG_ = self.___PRIVATE_empty_data_generator___
            self._KG_ = self.___PRIVATE_constant_cache_key_generator___

        else:
            try:
                _ = data_generator % 1
            except TypeError: # check `data_generator` is not a number!
                pass

            if cache_key_generator == 'no_cache':  # do not cache
                # use this when nothing is the same over all elements, we make the data whenever it is called.
                self._DG_ = data_generator
                self._KG_ = self.___PRIVATE_no_cache_key_generator___

            else:

                #If it's not no_cache, it must be a method (cannot be a function) have 1 input (#element) except self.
                if isinstance(cache_key_generator, types.MethodType):
                    # noinspection PyUnresolvedReferences
                    assert cache_key_generator.__code__.co_argcount == 2
                else:
                    pass

                self._DG_ = data_generator
                self._KG_ = cache_key_generator


        self._gathering_matrix_ = None
        self.RESET_cache()
        self.___CT___ = '>CT<'
        self.___NC___ = '>NC<'
        self.___IS_NC___ = False
        self.___IS_CT___ = False
        self.___CHECK_repeat_CT___ = True
        self.___repeat_CK___ = ''
        self._con_shape_ = con_shape
        self._customize_ = SpaVec_Customize(self)
        self._shape_ = None
        self._freeze_self_()

    def RESET_cache(self):
        self._cache_ = dict()

    def ___PRIVATE_empty_data_generator___(self, i):
        assert i in self.elements
        return spspa.csc_matrix((self.___EMPTY_LENGTH___,1))

    def ___PRIVATE_constant_cache_key_generator___(self, i):
        assert i in self.elements
        return self.___CT___

    def ___PRIVATE_no_cache_key_generator___(self, i):
        assert i in self.elements
        return self.___NC___

    @property
    def gathering_matrix(self):
        return self._gathering_matrix_

    @gathering_matrix.setter
    def gathering_matrix(self, gathering_matrix):
        if gathering_matrix.__class__.__name__ == 'Chain_Gathering_Matrix':
            pass
        else:
            if not isinstance(gathering_matrix, (list, tuple)):
                gathering_matrix = [gathering_matrix,]
            cgm0 = list()
            for _ in gathering_matrix:
                if _.__class__.__name__ == 'Gathering_Matrix':
                    cgm0.append(_)
                else:
                    cgm0.append(_.numbering.gathering)
            gathering_matrix = Chain_Gathering_Matrix(cgm0)
        self._gathering_matrix_ = gathering_matrix

        if self._shape_ is not None:
            assert self._gathering_matrix_.shape + (1,) == self._shape_

    @property
    def assembled(self):
        """
        Return an assembled global vector.

        :return:
        :rtype: GlobalVector
        """
        assert self.gathering_matrix is not None, "I have no gathering matrix"
        GI = self.gathering_matrix
        DEP = GI.GLOBAL_num_dofs
        ROW = list()
        DAT = list()
        for i in self:
            Vi = self[i]
            indices = Vi.indices
            data = Vi.data
            ROW.extend(GI[i][indices])
            DAT.extend(data)
        b = GlobalVector(spspa.csc_matrix((DAT, ROW, [0, len(ROW)]), shape=(DEP, 1)))

        assert b.shape == (GI.GLOBAL_num_dofs, 1)

        return b

    @property
    def elements(self):
        return self._elements_

    def __iter__(self):
        for i in self.elements:
            yield i

    def __getitem__(self, item):
        # ck = self._KG_(item)
        # if ck == self.___CT___:
        #     pass
        # if self.___NC___ in ck:
        #     RETURN = self._DG_(item)
        # else:
        #     # Current, we only no_cache vector (we do not really have situation that needs cache a vector.)
        #     raise NotImplementedError()
        #

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

        RETURN = self.customize.___PRIVATE_do_execute_customization___(RETURN, item)

        # here we do not check the shape of RETURN because we could use this method to compute shape, so if we do, we may end up in dead loop.

        return RETURN

    def __contains__(self, item):
        return item in self.elements

    def __len__(self):
        return len(self.elements)

    @property
    def con_shape(self):
        return self._con_shape_

    @property
    def shape(self):
        """The local shape : == (len(self),) + np.shape(self[i]). So the
        first value refer to how many local mesh elements. The second
        value refers to how many entries in the local vector. The third
        one must be 1.

        :return: A tuple of 3 integers. The third integer must be 1
        """
        if self._shape_ is not None: return self._shape_

        if self.gathering_matrix is not None:
            self._shape_ = self.gathering_matrix.shape + (1,)
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
                assert _s_[1] == 1, f"EWC vector must of shape (.,1)."
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
        """We use sub-methods of this properties to add customization. These customization will be
        executed when I am called."""
        return self._customize_

    def __sub__(self, other):
        """self - other"""
        assert other.__class__.__name__ == 'EWC_ColumnVector'
        assert self._elements_._mesh_ == other._elements_._mesh_
        DKC = ___CV_SUB___(self, other)
        return EWC_ColumnVector(self._elements_, DKC, DKC.__KG_call__)

    def __add__(self, other):
        """self + other"""
        assert other.__class__.__name__ == 'EWC_ColumnVector'
        assert self._elements_._mesh_ == other._elements_._mesh_
        DKC = ___CV_ADD___(self, other)
        return EWC_ColumnVector(self._elements_, DKC, DKC.__KG_call__)

    def __neg__(self):
        """- EWC_ColumnVector"""
        data_generator = ___CV_NEG___(self)
        RETURN = EWC_ColumnVector(self._elements_, data_generator, self._KG_)
        if self.gathering_matrix is not None:
            RETURN.gathering_matrix = self.gathering_matrix
        return RETURN


class SpaVec_Customize(FrozenOnly):
    def __init__(self, spa_vec):
        self._spa_vec_ = spa_vec
        self.___customizations___ = dict()
        self._freeze_self_()

    @property
    def _customizations_(self):
        """For example
        self.___customizations___ = {
                ...,
                5: [('cletr', 16), ...],  # The first customization for #5 element output is "clear local entry #16".
                    Customization are executed in a positive sequence.
                ...
        }
        """
        return self.___customizations___

    def ___PRIVATE_do_execute_customization___(self, RETURN, e):
        """Execute all added customization.

        :param RETURN: The RETURN to be customized.
        :param e: for #i element.
        :return:
        """
        if e not in self._customizations_:

            pass

        else:
            assert spspa.isspmatrix_csc(RETURN), "We need to start with a (1d) csc matrix."

            CUSi = self._customizations_[e] # customizations for the output of #e element.
            for cus in CUSi:

                key, factors = cus

                # clear a local row of the EWC-sparse-matrix ----------------- BELOW -----------------------------------
                if key == 'cletr': # Clear Local EnTRy #factors
                    # `factors` will be an int
                    if factors.__class__.__name__ in ("int", "int32", "int64"):
                        pass
                    elif isinstance(factors, float):
                        assert factors % 1 == 0, f"factors={factors} wrong, when `cletr`, factors must be an int."
                    else:
                        raise Exception(f"factors={factors} wrong, when `cletr`, factors must be an int.")

                    if not spspa.isspmatrix_lil(RETURN): RETURN = RETURN.tolil()
                    RETURN[factors, 0] = 0
                elif key == 'slest': # Set Local EntrieS To
                    if not spspa.isspmatrix_lil( RETURN): RETURN = RETURN.tolil()
                    RETURN[factors[0], 0] = factors[1]

                # Not Implemented --------------------------------------------- BELOW ----------------------------------
                else:
                    raise NotImplementedError(f"Can not handle customization key={key}.")
                #============================================================== ABOVE ==================================

        if not spspa.isspmatrix_csc(RETURN): RETURN = RETURN.tocsc()
        assert RETURN.shape[1] == 1, f"A vector must to csc_matrix of shape (x, 1), now its shape is {RETURN.shape}."
        return RETURN

    @accepts('self', (int, float, 'int32', 'int64'),
             ('PartialDofs', 'PartialCochain'), (int, float))
    def set_constant_entries_according_to_CSCG_partial_dofs(
            self, i, pd, constant, interpreted_as='local_dofs'):
        """ We do:
        block[i][row_pds.dofs, 1] = constant
        locally. Note it is "locally", that means, after assembly, we
        may have values which are multiple times of the constant. This is okay because same thing
        will happen in the left-hand-side EWC matrix.

        :param i:
        :param pd:
        :param constant:
        :param interpreted_as:
        :return:
        """

        if pd.__class__.__name__ == 'PartialCochain':
            pd = pd.dofs

        bsp = self._spa_vec_.con_shape
        assert bsp is not False and np.shape(bsp) == (1,), \
            "we can only use this to concatenate EWC vectors to keep the block safe. " \
            "if you want to apply it to a single EWC vectors, first wrap it in a list, then send it to concatenate"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare vector of block shape {bsp}."

        if interpreted_as == 'local_dofs':
            GM = self._spa_vec_.gathering_matrix
            LDR_row = GM.local_dofs_ranges[i]
            start_row = LDR_row.start
            dofs = pd.interpreted_as.local_dofs
            for e in dofs:
                local_dofs = np.array(dofs[e]) + start_row

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                self.___customizations___[e].append(('slest', (local_dofs, constant)))

        else:
            raise NotImplementedError(f"interpreted_as={interpreted_as} not implemented.")



    @accepts('self', (int, float, 'int32', 'int64'), ('PartialDofs', 'PartialCochain'))
    def set_entries_according_to_two_CSCG_partial_cochains(
        self, i, pd, pc=None, interpreted_as='local_dofs'):
        """We do:
        block[i][row_pds.dofs, 1] = col_pds.cochain
        locally. Note it is "locally", that means, after assembly, we
        may have values which are multiple times of the value. This is okay because same thing
        will happen in the left-hand-side EWC matrix.

        :param i:
        :param pd:
        :param pc:
        :param interpreted_as:
        :return:
        """
        if pc is None:
            assert pd.__class__.__name__ == 'PartialCochain', \
                "I need a PartialCochain when pc is None."
            pc = pd

        if pd.__class__.__name__ == 'PartialCochain': pd = pd.dofs
        assert pd.__class__.__name__ == 'PartialDofs', f"pd must be a PartialDofs."
        assert pc.__class__.__name__ == 'PartialCochain', f"pc must be a PartialCochain."

        bsp = self._spa_vec_.con_shape
        assert bsp is not False and np.shape(bsp) == (1,), \
            "we can only use this to concatenate EWC vectors to keep the block safe. " \
            "if you want to apply it to a single EWC vectors, first wrap it in a list, then send it to concatenate"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare vector of block shape {bsp}."


        if interpreted_as == 'local_dofs':
            GM = self._spa_vec_.gathering_matrix
            LDR_row = GM.local_dofs_ranges[i]
            start_row = LDR_row.start
            dofs = pd.interpreted_as.local_dofs
            cochain = pc.cochain
            assert len(dofs) == len(cochain), "pc, pd must cover same amount of mesh elements."
            for e in dofs:
                local_dofs = np.array(dofs[e]) + start_row
                local_cochain = cochain[e]
                if e not in self.___customizations___:
                    self.___customizations___[e] = list()
                # Set Local Entries To
                self.___customizations___[e].append(('slest', (local_dofs, local_cochain)))

        else:
            raise NotImplementedError(f"interpreted_as={interpreted_as} not implemented.")




class ___CV_SUB___(FrozenOnly):
    def __init__(self, v1, v2):
        self._v1_ = v1
        self._v2_ = v2
        self._freeze_self_()

    def __call__(self, i):
        """"""
        return self._v1_[i] - self._v2_[i]

    def __KG_call__(self, i):
        """"""
        return self._v1_._KG_(i) + self._v1_._KG_(i)

class ___CV_ADD___(FrozenOnly):
    def __init__(self, v1, v2):
        self._v1_ = v1
        self._v2_ = v2
        self._freeze_self_()

    def __call__(self, i):
        """"""
        return self._v1_[i] + self._v2_[i]

    def __KG_call__(self, i):
        """"""
        return self._v1_._KG_(i) + self._v1_._KG_(i)

class ___CV_NEG___(FrozenOnly):
    def __init__(self, V):
        self._V_ = V
        self._freeze_self_()

    def __call__(self, i):
        return -self._V_[i]











def ___concatenate_EWC_sparse_vectors___(vectors):
    """"""
    assert isinstance(vectors, (list, tuple)) and len(vectors) > 0, "please put (more than 0) vectors in list or tuple."
    elements = None
    for i, v in enumerate(vectors):
        assert v.__class__.__name__ == 'EWC_ColumnVector', f"Cannot handle vectors={vectors}."
        assert v.customize._customizations_ == dict(), \
            f"vectors[{i}] is customized, can not used for concatenate."

        if elements is None:
            elements = v.elements
        else:
            assert elements == v.elements

    DG = ___concatenate_HELPER_DataGenerator___(vectors)
    KG = ___concatenate_HELPER_KeyGenerator___(vectors)

    EWC = EWC_ColumnVector(elements, DG, KG, con_shape=(len(vectors),))

    CGMs = list()
    for v in vectors:
        CGMs.append(v.gathering_matrix)

    if None not in CGMs:
        EWC.gathering_matrix = Chain_Gathering_Matrix(CGMs)

    return EWC

class ___concatenate_HELPER_DataGenerator___:
    """"""
    def __init__(self, vectors):
        self.vectors = vectors
        self.I = len(vectors)

    def __call__(self, i):
        """"""
        output = [None for _ in range(self.I)]
        for j, Vj in enumerate(self.vectors):
            output[j] = Vj[i]
        return spspa.vstack(output, format='csc')

class ___concatenate_HELPER_KeyGenerator___:
    """"""
    def __init__(self, vectors):
        self.vectors = vectors

    def __call__(self, i):
        """"""
        key_out = ''
        for v in self.vectors:
            key_out += v._KG_(i)
        return key_out














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
        self.RESET_cache()
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


    def RESET_cache(self):
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
        """We use sub-methods of this properties to add customization. These customization will be
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



class SpaMat_Customize(FrozenOnly):
    """Store all method to customize the EWC sparse matrix."""
    def __init__(self, spa_mat):
        self._spa_mat_ = spa_mat
        self.___customizations___ = dict()
        self._freeze_self_()

    @property
    def _customizations_(self):
        """For example
        self.___customizations___ = {
                ...,
                5: [('clr', 16), ...],  # The first customization for #5 element output is "clear local row #16".
                    Customization are executed in a positive sequence.
                6: [('salv', ([5, 6], 151.1256221)), ...],  # Set A Local Value whose local indices are (5,6) in element 6 to
                    151.1256221.
                ...
        }
        """
        return self.___customizations___

    def ___PRIVATE_do_execute_customization___(self, RETURN, e):
        """Execute all added customization.

        :param RETURN: The RETURN to be customized.
        :param e: for #e element.
        :return:
        """

        # print('CUSTOMIZATION')

        if e not in self._customizations_:

            pass

        else:
            assert spspa.isspmatrix_csc(RETURN) or spspa.isspmatrix_csr(RETURN)

            CUSi = self._customizations_[e] # customizations for the output of #e element.



            for cus in CUSi:

                key, factors = cus # for example cus=('clr', 5) or ('salv', ([1,1],5))

                # clear a local row of the EWC-sparse-matrix ----------------- BELOW -----------------------------------
                if key == 'clr': # Clear Local Row
                    # `factors` will be an int
                    assert factors % 1 == 0, f"factors={factors} wrong, when `clr`, factors must be an int."
                    if not spspa.isspmatrix_lil(RETURN): RETURN = RETURN.tolil()
                    RETURN[factors, :] = 0

                # clear local rows for element #e --------------------- BELOW ---------------------------
                elif key == 'clrs': # Clear Local Row
                    #factores are rows
                    if not spspa.isspmatrix_lil(RETURN): RETURN = RETURN.tolil()
                    RETURN[factors, :] = 0

                # Set A Local Value of local indices (i, j) to v --------------------- BELOW ---------------------------
                elif key == 'salv':
                    ij, v = factors
                    i, j = ij
                    if not spspa.isspmatrix_lil(RETURN): RETURN = RETURN.tolil()
                    RETURN[i, j] = v

                # identify Local Rows --------------------- BELOW ---------------------------
                elif key == 'ilrs':
                    if not spspa.isspmatrix_lil( RETURN): RETURN = RETURN.tolil()
                    RETURN[factors, :] = 0
                    RETURN[factors, factors] = 1

                # identify Local Rows At Columns --------------------- BELOW ---------------------------
                elif key == 'ilrsac':
                    # key, factors = ('ilrsac', ([0,1,2,3], [10, 11, 12 ,14]))
                    if not spspa.isspmatrix_lil( RETURN): RETURN = RETURN.tolil()
                    RETURN[factors[0], :] = 0
                    RETURN[factors[0], factors[1]] = 1
                    # print(e, key, factors)
                    # print(len(factors[0]))


                # Not Implemented --------------------------------------------- BELOW ----------------------------------
                else:
                    raise NotImplementedError(f"Can not handle customization key={key}.")
                #============================================================== ABOVE ==================================


        if not (spspa.isspmatrix_csc(RETURN) or spspa.isspmatrix_csr(RETURN)): RETURN = RETURN.tocsr()


        return RETURN


    @accepts('self', (int, float, 'int32', 'int64'))
    def clear_global_row(self, r):
        """Make the #r row of the global matrix (assemble self to get the global matrix) to be all zero.
        To achieve this, we only adjust the local output of the call function.

        We need to first find which elements contribute to global row #r, and then find the "local rows".
        Afterwards, we can apply customization key 'clr' with customization factor "local row" to elements.

        :param r:
        :return:
        """
        assert r % 1 == 0, f"r={r}({r.__class__.__name__}) is wrong."
        if isinstance(r, float): r = int(r)

        rCGM = self._spa_mat_.gathering_matrices[0]
        assert rCGM is not None, "I have no row gathering_matrix!"
        OUTPUT = rCGM.FIND.elements_and_local_indices_of_dof(r)
        if OUTPUT is None:
            pass
        else:
            elements, local_indices = OUTPUT
            for e, i in zip(elements, local_indices):

                assert e in self._spa_mat_, f"element {e} is not in this core!"

                if e not in self._customizations_:
                    self._customizations_[e] = list()

                self._customizations_[e].append(('clr', int(i)))

    def clear_global_rows(self, rs):
        """

        :param rs:
        :return:
        """
        raise NotImplementedError()


    @accepts('self', (int, float, 'int32', 'int64'), (int, float, 'int32', 'int64'), (int, float, 'int32', 'int64'))
    def set_assembled_M_ij_to(self, i, j, v):
        """Let M be the assembled matrix, we set M[i,j] = v.

        :param i:
        :param j:
        :param v:
        :return:
        """
        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) is wrong."
        if isinstance(i, float): i = int(i)
        assert j % 1 == 0, f"j={j}({j.__class__.__name__}) is wrong."
        if isinstance(j, float): j = int(j)

        rCGM, cCGM = self._spa_mat_.gathering_matrices
        assert rCGM is not None, "I have no row gathering_matrix!"
        assert cCGM is not None, "I have no col gathering_matrix!"

        rO = rCGM.FIND.elements_and_local_indices_of_dof(i)
        cO = cCGM.FIND.elements_and_local_indices_of_dof(j)

        if rO is not None and cO is not None:
            rE, rI = rO
            cE, cI = cO

            TBG = list()

            LOCAL_ij = list()

            for i, re in enumerate(rE):
                for j, ce in enumerate(cE):
                    if ce == re:
                        TBG.append(re)
                        LOCAL_ij.append((rI[i], cI[j]))
                        break
                    else:
                        pass

            if TBG == list():
                TBG = None

        else:
            TBG = None # to be gathered to master

        ALL_TBG = cOmm.gather(TBG, root=mAster_rank)
        if rAnk == mAster_rank:
            ELEMENTS_HAVE_THE_TARGET = list()
            for tbg in ALL_TBG:
                if tbg is not None:
                    ELEMENTS_HAVE_THE_TARGET.extend(tbg)
                else:
                    pass
            ELEMENTS_HAVE_THE_TARGET = list(set(ELEMENTS_HAVE_THE_TARGET))

            assert len(ELEMENTS_HAVE_THE_TARGET) > 0, "We are trying to set value to an empty place at which the " \
                                                      "assembling will allocate no value."


            the_element = ELEMENTS_HAVE_THE_TARGET[0] # we must find an element by here.
        else:
            the_element = None
        the_element = cOmm.bcast(the_element, root=mAster_rank)


        if TBG is not None:
            if the_element in TBG:
                I_am_in = rAnk
            else:
                I_am_in = sIze
        else:
            I_am_in = sIze

        who_are_in = cOmm.gather(I_am_in, root=mAster_rank)
        if rAnk == mAster_rank:
            use_who = min(who_are_in)
            assert use_who != sIze, "Something is wrong, we need to find a correct core."
        else:
            use_who = None
        use_who = cOmm.bcast(use_who, root=mAster_rank)

        if use_who == rAnk:
            assert the_element in TBG, "Must be like this!"
            for k, e in enumerate(TBG):
                # noinspection PyUnboundLocalVariable
                i, j = LOCAL_ij[k]
                if e not in self.___customizations___:
                    self.___customizations___[e] = list()
                if e != the_element:
                    self.___customizations___[e].append(('salv', ([i, j], 0)))
                else:
                    self.___customizations___[e].append(('salv', ([i, j], v)))



        else: # These cores has no business with the setting the value, they only need to make the value to be zero.
            if TBG is None: # these cores have no business at all.
                pass
            else:
                for k, e in enumerate(TBG):
                    # noinspection PyUnboundLocalVariable
                    i, j = LOCAL_ij[k]
                    if e not in self.___customizations___:
                        self.___customizations___[e] = list()
                    self.___customizations___[e].append(('salv', ([i, j], 0)))


    @accepts('self', (int, float, 'int32', 'int64'))
    def identify_global_row(self, r):
        """Let M be the assembled matrix. We want to make M[r,:] = 0 except M[r,r] = 1.

        :param r:
        :return:
        """
        assert r % 1 == 0, f"r={r}({r.__class__.__name__}) is wrong."
        if isinstance(r, float): r = int(r)
        self.clear_global_row(r)
        self.set_assembled_M_ij_to(r, r, 1)


    def identify_global_rows(self, rs):
        """

        :param rs:
        :return:
        """
        raise NotImplementedError()


    @accepts('self', (int, float, 'int32', 'int64'), ('PartialDofs', 'PartialCochain'))
    def identify_global_rows_according_to_CSCG_partial_dofs(self, i, pds, interpreted_as='local_dofs'):
        """We locally make M[r,:] = 0 except M[r,r] = 1 according to a CSCG PartialDofs
        instance in block[i][:]

        This will ask the local and global matrix to be square: But we
        do not check it.

        IMPORTANT: as this is done only locally, after assembly, M[r,r] may be greater than 1.
            for example, for a corner point, M[r,r] can be 8 if it is internal. This, for example,
            will not introduce problem for imposing bc if in the right-hand-side vector, we
            also give a 8 times greater value.

        :param pds: we will get dofs form this.
        :param i: We locally make block[i][:][dofs, :] = 0 except block[i][i][dofs, dofs] = 1
        :param interpreted_as:
        :return:
        """
        if pds.__class__.__name__ == 'PartialCochain':
            pds = pds.dofs

        assert self._spa_mat_.elements._mesh_ == pds._mesh_, \
            "PartialDofs elements != EWC elements."


        bsp = self._spa_mat_.bmat_shape

        assert bsp is not False and np.shape(bsp) == (2,), \
            "we can only use this to bmat EWC matrices to keep the block safe. " \
            "if you want to apply it to a single EWC matrix, first wrap it in a list, then send it to bmat"

        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare matrix of block shape ({bsp[0]},.)."

        assert self._spa_mat_._gathering_matrices_0_ is not None, \
            f"We at least need row gathering matrix!"


        if interpreted_as == 'local_dofs':

            GM = self._spa_mat_.gathering_matrices[0]
            LDR = GM.local_dofs_ranges[i]
            start = LDR.start

            LDF = pds.interpreted_as.local_dofs

            for e in LDF: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                local_dofs = np.array(LDF[e]) + start

                self.___customizations___[e].append(('ilrs', local_dofs))

        else:
            raise Exception(f"Cannot identify global rows through "
                            f"SCG_partial_dofs interpreted "
                            f"as <{interpreted_as}>.")



    @accepts('self', (int, float, 'int32', 'int64'), (int, float, 'int32', 'int64'),
             ('PartialDofs', 'PartialCochain'),('PartialDofs', 'PartialCochain'))
    def off_diagonally_identify_rows_according_to_two_CSCG_partial_dofs(
        self, i, j, row_pds, col_pds, interpreted_as='local_dofs'):
        """We will identify off-diagonal block, block[i][j], and set
        block[i][:][row_pds.dofs, :] = 0
        block[i][j][row_pds.dofs, col_pds.dofs] = 1
        locally. Note it is "locally", that means, after assembly, we
        may have values greater than 1. This is okay because same thing
        will happen in the right-hand-side EWC vector.

        :param i:
        :param j:
        :param row_pds:
        :param col_pds:
        :param interpreted_as:
        :return:
        """
        if row_pds.__class__.__name__ == 'PartialCochain':
            row_pds = row_pds.dofs
        if col_pds.__class__.__name__ == 'PartialCochain':
            col_pds = col_pds.dofs

        assert self._spa_mat_.elements._mesh_ == row_pds._mesh_, \
            "row PartialDofs elements != EWC elements."
        assert self._spa_mat_.elements._mesh_ == col_pds._mesh_, \
            "col PartialDofs elements != EWC elements."


        assert self._spa_mat_._gathering_matrices_0_ is not None, \
            f"We need row gathering matrix!"
        assert self._spa_mat_._gathering_matrices_1_ is not None, \
            f"We need col gathering matrix!"

        bsp = self._spa_mat_.bmat_shape
        assert bsp is not False and np.shape(bsp) == (2,), \
            "we can only use this to bmat EWC matrices to keep the block safe. " \
            "if you want to apply it to a single EWC matrix, first wrap it in a list, then send it to bmat"


        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        if isinstance(i, float): i = int(i)
        assert 0 <= i < bsp[0], f"i={i} is beyond the shape of this " \
                                f"EWC spare matrix of block shape ({bsp[0]},.)."

        assert j % 1 == 0, f"j={j}({j.__class__.__name__}) cannot be an index!"
        if isinstance(j, float): j = int(j)
        assert 0 <= j < bsp[1], f"j={j} is beyond the shape of this " \
                                f"EWC spare matrix of block shape ({bsp[1]},.)."

        if interpreted_as == 'local_dofs':

            GM = self._spa_mat_.gathering_matrices
            LDR_row = GM[0].local_dofs_ranges[i]
            LDR_col = GM[1].local_dofs_ranges[j]
            start_row = LDR_row.start
            start_col = LDR_col.start

            LDF_row = row_pds.interpreted_as.local_dofs
            LDF_col = col_pds.interpreted_as.local_dofs

            # print(GM[0].local_dofs_ranges, start_row, start_col)

            for e in LDF_row: # go through all locally involved mesh element numbers
                assert e in self._spa_mat_ and e in LDF_col

                if e not in self.___customizations___:
                    self.___customizations___[e] = list()

                row_local_dofs = np.array(LDF_row[e]) + start_row
                col_local_dofs = np.array(LDF_col[e]) + start_col

                self.___customizations___[e].append(
                    ('ilrsac', (row_local_dofs, col_local_dofs)))

                # print(row_local_dofs, col_local_dofs)

        else:
            raise Exception(f"Cannot off-diagonally identify global rows through "
                            f"SCG_partial_dofs interpreted "
                            f"as <{interpreted_as}>.")





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


















def ___bmat_EWC_sparse_matrices___(blocks):
    """A function to do bmat of EWC_SparseMatrix.

    :param blocks:
    :return:
    """
    assert isinstance(blocks, (list, tuple)), "please put blocks in list or tuple."

    I = len(blocks)
    J = None
    for i, bR in enumerate(blocks):
        assert isinstance(bR, (list, tuple)), "please put blocks in list or tuple."
        if J is None:
            J = len(bR)
        else:
            assert J == len(bR)

    elements = None

    RGM = [[None for _ in range(J)] for _ in range(I)]
    CGM = [[None for _ in range(J)] for _ in range(I)]

    for i, Bi in enumerate(blocks):
        for j, Bij in enumerate(Bi):
            if Bij is None:
                RGM[i][j] = None
                CGM[i][j] = None
            else:
                assert Bij.customize._customizations_ == dict(), \
                    f"customized block[{i}][{j}] cannot be used for bmat."
                # this is because of the cache function. the customizations of block will not be renewed in cache.

                assert Bij.__class__.__name__ == 'EWC_SparseMatrix', \
                    f"I only handle EWC_SparseMatrix, but block[{i}][{j}] now is {Bij.__class__.__name__}."

                if elements is None:
                    elements = Bij.elements
                else:
                    assert Bij.elements == elements

                RGM[i][j], CGM[i][j] = Bij.gathering_matrices

    assert elements is not None, f"blocks of only None?"
    DG = ___BMAT_HELPER_DataGenerator___(blocks)
    KG = ___BMAT_HELPER_KeyGenerator___(blocks)
    EWC = EWC_SparseMatrix(elements, DG, KG, bmat_shape=(I, J))

    # Now we have a look at if we can get the gathering matrices for the bmat result.
    rgm = [None for _ in range(I)]
    for i in range(I):
        for j in range(J):
            Rij = RGM[i][j]
            if Rij is None:
                pass
            else:
                # noinspection PyUnresolvedReferences
                assert str(Rij.__class__) == "<class 'TOOLS.linear_algebra.gathering.Chain_Gathering_Matrix'>"
                if rgm[i] is None:
                    rgm[i] = Rij
                else:
                    assert rgm[i] == Rij, f"bmat wrong blocks. Gathering matrices for row[{i}] dis-match."

    cgm = [None for _ in range(J)]
    for j in range(J):
        for i in range(I):
            Cij = CGM[i][j]
            if Cij is None:
                pass
            else:
                # noinspection PyUnresolvedReferences
                assert str(Cij.__class__) == "<class 'TOOLS.linear_algebra.gathering.Chain_Gathering_Matrix'>"
                if cgm[j] is None:
                    cgm[j] = Cij
                else:
                    assert cgm[j] == Cij, f"bmat wrong blocks. Gathering matrices for col[{j}] dis-match."

    doR = False if None in rgm else True
    doC = False if None in cgm else True

    if doR and doC:
        R_CGM = Chain_Gathering_Matrix(rgm)
        C_CGM = Chain_Gathering_Matrix(cgm)
        EWC.gathering_matrices = (R_CGM, C_CGM)

    if 0 in EWC: _ = EWC[0]
    # do a test to check if bmat is fine. this is OKAY even if we will apply
    # some customization later, the cache is done before the customization

    return EWC

class ___BMAT_HELPER_DataGenerator___:
    """"""
    def __init__(self, blocks):
        self.blocks = blocks
        self.I = len(blocks)
        self.J = len(blocks[0])

    def __call__(self, item):
        """"""
        output = [[None for _ in range(self.J)] for _ in range(self.I)]
        for i, Bi in enumerate(self.blocks):
            for j, Bij in enumerate(Bi):
                if Bij is not None:
                    output[i][j] = Bij[item]
        return spspa.bmat(output, format='csr')

class ___BMAT_HELPER_KeyGenerator___:
    """"""
    def __init__(self, blocks):
        self.valid_block = list()
        for i, Bi in enumerate(blocks):
            for j, Bij in enumerate(Bi):
                if Bij is not None:
                    self.valid_block.append(Bij)

    def __call__(self, item):
        key_out = ''
        for vb in self.valid_block:
            key_out += vb._KG_(item)
        return key_out