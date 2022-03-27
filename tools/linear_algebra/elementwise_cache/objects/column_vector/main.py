


import types
from screws.freeze.main import FrozenOnly
from scipy import sparse as spspa
from root.config.main import *
from tools.linear_algebra.gathering.chain_matrix.main import Chain_Gathering_Matrix
from tools.linear_algebra.data_structures.global_matrix.main import GlobalVector

from tools.linear_algebra.elementwise_cache.objects.column_vector.helpers.add import ___CV_ADD___
from tools.linear_algebra.elementwise_cache.objects.column_vector.helpers.neg import ___CV_NEG___
from tools.linear_algebra.elementwise_cache.objects.column_vector.helpers.sub import ___CV_SUB___

from tools.linear_algebra.elementwise_cache.objects.column_vector.customize import SpaVec_Customize
from tools.linear_algebra.elementwise_cache.objects.column_vector.assembler import EWC_ColumnVector_Assembler

class EWC_ColumnVector(FrozenOnly):
    """
    Element-wise cached column vector (csc_matrix of shape (n,1)).

    example::

        - EWC_ColumnVector(mesh, 5) : generate an empty local vector of shape (5,1) for all elements.
        - EWC_ColumnVector(mesh, form): make an empty local vector of shape (form.num.basis, 1).
            So, equivalent to EWC_ColumnVector(mesh, form.num.basis)

    :param mesh_elements:
    :param data_generator:
    :param con_shape:
    :type con_shape: bool, tuple, list
    """
    def __init__(self, mesh_elements, data_generator, cache_key_generator=None, con_shape=False):
        """

        :param mesh_elements:
        :param data_generator: {}
        :param cache_key_generator:
        :param con_shape:
        """
        if mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh_Elements', '_2dCSCG_Mesh_Elements'):
            self._elements_ = mesh_elements
        elif mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'):
            self._elements_ = mesh_elements.elements
        else:
            raise Exception()

        # --------- parse data type ---------------------------------------------------------------
        DATA_TYPE = None

        if data_generator.__class__.__name__ in ('int', 'float', 'int32', 'int64'):
            DATA_TYPE = 'EMPTY'
        elif hasattr(data_generator, 'standard_properties') and \
            'CSCG_form' in data_generator.standard_properties.tags:
            # data_generator is a CSCG form.
            DATA_TYPE = 'EMPTY'
            data_generator = data_generator.num.basis
            # equivalent to EWC_ColumnVector(mesh, form.num.basis)
        else: # default data type
            pass

        #---------- parse default cache_key_generator ---------------------------------------------
        if cache_key_generator is None:

            if DATA_TYPE == 'EMPTY': # make empty sparse vector.

                # noinspection PyTypeChecker
                if data_generator % 1 == 0:

                    self.___EMPTY_LENGTH___ = data_generator

                else:
                    raise Exception(f"empty sparse vector of local length={data_generator} is wrong.")

            else:
                cache_key_generator = 'no_cache' # do NOT CHANGE THIS DEFAULT SETTING!!!

        else: # have given a generator, use it.
            pass

        #--------- get DG and KG -------------------------------------------------------------------

        if DATA_TYPE == 'EMPTY':

            self._DG_ = self.___PRIVATE_empty_data_generator___
            self._KG_ = self.___PRIVATE_constant_cache_key_generator___

        elif DATA_TYPE is None:

            assert callable(data_generator), f"data_generator={data_generator} is not callable!"

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

        else:
            raise Exception(f"DATA_TYPE = {DATA_TYPE} wrong!")

        #--------------------------------------- DONE ----------------------------------------

        self._gathering_matrix_ = None
        self.___PRIVATE_reset_cache___()
        self.___CT___ = '>CT<'
        self.___NC___ = '>NC<'
        self.___IS_NC___ = False
        self.___IS_CT___ = False
        self.____CT_DG____ = None
        self.___CHECK_repeat_CT___ = True
        self.___repeat_CK___ = ''
        self._con_shape_ = con_shape
        self._customize_ = SpaVec_Customize(self)
        self._shape_ = None
        self._assembler_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._cache_ = dict()

    def ___PRIVATE_empty_data_generator___(self, i):
        assert i in self.elements
        return spspa.csc_matrix((self.___EMPTY_LENGTH___, 1))

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
    def assembler(self):
        if self._assembler_ is None:
            self._assembler_ = EWC_ColumnVector_Assembler(self)
        return self._assembler_

    @property
    def assembled(self):
        """
        Return an assembled global vector.

        :return:
        :rtype: GlobalVector
        """
        return self.assembler()

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
        """We use sub-methods of these properties to add customization. These customizations will be
        executed when I am called.
        """
        return self._customize_

    def __sub__(self, other):
        """self - other"""
        assert other.__class__.__name__ == 'EWC_ColumnVector'
        assert self._elements_._mesh_ == other._elements_._mesh_
        DKC = ___CV_SUB___(self, other)
        # noinspection PyTypeChecker
        return EWC_ColumnVector(self._elements_, DKC, DKC.__KG_call__)

    def __add__(self, other):
        """self + other"""
        assert other.__class__.__name__ == 'EWC_ColumnVector'
        assert self._elements_._mesh_ == other._elements_._mesh_
        DKC = ___CV_ADD___(self, other)
        # noinspection PyTypeChecker
        return EWC_ColumnVector(self._elements_, DKC, DKC.__KG_call__)

    def __neg__(self):
        """- EWC_ColumnVector"""
        data_generator = ___CV_NEG___(self)
        # noinspection PyTypeChecker
        RETURN = EWC_ColumnVector(self._elements_, data_generator, self._KG_)
        if self.gathering_matrix is not None:
            RETURN.gathering_matrix = self.gathering_matrix
        return RETURN