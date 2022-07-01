# -*- coding: utf-8 -*-


import types
from screws.freeze.main import FrozenOnly
from scipy import sparse as spspa
from tools.linear_algebra.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix
from tools.linear_algebra.gathering.irregular.ir_chain_matrix.main import iR_Chain_Gathering_Matrix
from tools.linear_algebra.data_structures.global_matrix.main import GlobalVector

from tools.linear_algebra.elementwise_cache.objects.column_vector.helpers.add import ___CV_ADD___
from tools.linear_algebra.elementwise_cache.objects.column_vector.helpers.neg import ___CV_NEG___
from tools.linear_algebra.elementwise_cache.objects.column_vector.helpers.sub import ___CV_SUB___

from tools.linear_algebra.elementwise_cache.objects.column_vector.customize import SpaVec_Customize
from tools.linear_algebra.elementwise_cache.objects.column_vector.assembler import EWC_ColumnVector_Assembler
from tools.linear_algebra.elementwise_cache.objects.column_vector.adjust import SpaVec_Adjust
from tools.linear_algebra.elementwise_cache.objects.column_vector.blocks.main import EWC_ColVec_Blocks

class EWC_ColumnVector(FrozenOnly):
    """
    Element-wise cached column vector (csc_matrix of shape (n,1)).

    example::

        - EWC_ColumnVector(mesh, 5) : generate an empty local vector of shape (5,1) for all elements.

        - EWC_ColumnVector(mesh, form): make an empty local vector of shape (form.num.basis, 1).
            So, equivalent to EWC_ColumnVector(mesh, form.num.basis).

            Note that this firm can be an AD-form.

        - EWC_ColumnVector(form) same as EWC_ColumnVector(form.mesh, form)

        - EWC_ColumnVector(mesh, dict): the data are in a dict already.


    :param mesh_elements:
    :param data_generator:
    :param con_shape:
    :type con_shape: bool, tuple, list
    """
    def __init__(self, mesh_elements, data_generator=None, cache_key_generator=None, con_shape=False):
        """

        :param mesh_elements:
        :param data_generator:
        :param cache_key_generator:
        :param con_shape:
        """
        #----------- when only provide one input `mesh_elements` ------------------------------
        if data_generator is None and cache_key_generator is None and con_shape is False:

            #------- we have a cscg form -------------------------------------------------
            if hasattr(mesh_elements, 'standard_properties') and \
                'CSCG_form' in mesh_elements.standard_properties.tags:
                data_generator = mesh_elements
                mesh_elements = mesh_elements.mesh.elements
            elif hasattr(mesh_elements, 'prime') and \
                hasattr(mesh_elements.prime, 'standard_properties') and \
                'CSCG_form' in mesh_elements.prime.standard_properties.tags:
                data_generator = mesh_elements.prime
                mesh_elements = mesh_elements.mesh.elements

            #----------- we have a mpRfT form as the only input -------------------------------------
            elif hasattr(mesh_elements, 'standard_properties') and \
                'mpRfT_form' in mesh_elements.standard_properties.tags:
                self._mpRfT_num_basis_ = mesh_elements.num.basis
                data_generator = self.___Pr_mpRfT_empty_DG___
                mesh_elements = mesh_elements.mesh


        else:
            pass

        #------------------------ check ----------------------------------------------------------
        if mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh_Elements', '_2dCSCG_Mesh_Elements'):
            self._elements_ = mesh_elements
        elif mesh_elements.__class__.__name__ in ('_3dCSCG_Mesh', '_2dCSCG_Mesh'):
            self._elements_ = mesh_elements.elements
        elif mesh_elements.__class__.__name__ in ('mpRfT2_Mesh', 'mpRfT3_Mesh'):
            self._elements_ = mesh_elements.rcfc
        elif isinstance(mesh_elements, dict):
            self._elements_ = mesh_elements
        else:
            raise Exception(f"mesh_elements={mesh_elements} not understandable!")


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
        elif hasattr(data_generator, 'prime') and \
            hasattr(data_generator.prime, 'standard_properties') and \
            'CSCG_form' in data_generator.prime.standard_properties.tags:
            # data_generator is a CSCG AD-form.
            DATA_TYPE = 'EMPTY'
            data_generator = data_generator.prime.num.basis
            # equivalent to EWC_ColumnVector(mesh, ad_form.prime.num.basis)
        else: # default data type
            pass

        # we can accept a dictionary as a data generator, we will wrap it with a method -----------
        if isinstance(data_generator, dict):
            assert len(data_generator) == len(self._elements_), "dict key wrong."
            for _ in data_generator: assert _ in self._elements_, "dict key wrong."
            self.___dict_DG___ = data_generator
            data_generator = self.___PRIVATE_dict_2_method_data_generator___
            if cache_key_generator is None:
                cache_key_generator = 'no_cache'
                # the data are in the dict anyway, so do not need to be cached.
        else:
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
        self.___mpRfT_empty_cache___ = dict()
        self.___CHECK_repeat_CT___ = True
        self.___repeat_CK___ = ''
        self._con_shape_ = con_shape
        self._customize_ = SpaVec_Customize(self)
        self._assembler_ = None
        self._adjust_ = None
        self._blocks_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._cache_ = dict()


    def ___Pr_mpRfT_empty_DG___(self, rc_rp):
        """"""
        num_basis = self._mpRfT_num_basis_[rc_rp]

        if num_basis in self.___mpRfT_empty_cache___:
            return self.___mpRfT_empty_cache___[num_basis]
        else:
            EP = spspa.csc_matrix((num_basis, 1))
            self.___mpRfT_empty_cache___[num_basis] = EP
            return EP


    def ___PRIVATE_empty_data_generator___(self, i):
        assert i in self.elements
        return spspa.csc_matrix((self.___EMPTY_LENGTH___, 1))

    def ___PRIVATE_constant_cache_key_generator___(self, i):
        assert i in self.elements
        return self.___CT___

    def ___PRIVATE_no_cache_key_generator___(self, i):
        assert i in self.elements
        return self.___NC___

    def ___PRIVATE_dict_2_method_data_generator___(self, i):
        """When we get a dict as the data generator, we make wrap it with a method."""
        return self.___dict_DG___[i]

    @property
    def gathering_matrix(self):
        return self._gathering_matrix_

    @gathering_matrix.setter
    def gathering_matrix(self, gathering_matrix):
        """"""
        if gathering_matrix.__class__.__name__ in ('Chain_Gathering_Matrix', 'iR_Chain_Gathering_Matrix'):
            pass
        else:
            if not isinstance(gathering_matrix, (list, tuple)):
                gathering_matrix = [gathering_matrix,]

            cgm0 = list()
            for _ in gathering_matrix:
                if _.__class__.__name__ in ('Gathering_Matrix', 'iR_Gathering_Matrix,'):
                    cgm0.append(_)
                else:
                    if hasattr(_, '___IS_ADF___') and _.___IS_ADF___:
                        cgm0.append(_.prime.numbering.gathering)
                    else:
                        cgm0.append(_.numbering.gathering)

            if any([_.__class__.__name__ == 'iR_Gathering_Matrix' for _ in cgm0]):
                gathering_matrix = iR_Chain_Gathering_Matrix(cgm0)
            else:
                gathering_matrix = Chain_Gathering_Matrix(cgm0)
        self._gathering_matrix_ = gathering_matrix

    @property
    def adjust(self):
        if self._adjust_ is None:
            self._adjust_ = SpaVec_Adjust(self)
        return self._adjust_

    @property
    def blocks(self):
        if self._blocks_ is None:
            self._blocks_ = EWC_ColVec_Blocks(self)
        return self._blocks_

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
        # ck = self._KG_(i)
        # if ck == self.___CT___:
        #     pass
        # if self.___NC___ in ck:
        #     RETURN = self._DG_(i)
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
        DKC = ___CV_SUB___(self, other)
        # noinspection PyTypeChecker
        return EWC_ColumnVector(self._elements_, DKC, DKC.__KG_call__)

    def __add__(self, other):
        """self + other"""
        assert other.__class__.__name__ == 'EWC_ColumnVector'
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