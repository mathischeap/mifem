# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/7/18 22:19
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
import numpy as np

from tools.linearAlgebra.elementwiseCache.objects.multiDimMatrix.do.main import MDM_Do
from tools.linearAlgebra.elementwiseCache.objects.multiDimMatrix.IS import MDM_IS

from tools.linearAlgebra.elementwiseCache.objects.multiDimMatrix.helpers.__neg__ import ___NEG___
from tools.linearAlgebra.elementwiseCache.objects.multiDimMatrix.helpers.__mul__ import MDM_MUL

class MultiDimMatrix(FrozenOnly):
    def __init__(self, iterator, mdm_generator, correspondence, cache_key_generator):
        """"""
        self._iterator_ = iterator

        if isinstance(mdm_generator, dict):
            self._dict_mdm_data_ = mdm_generator
            self._DG_ = self.___dict_DG___
            assert cache_key_generator == 'no_cache'
            self._KG_ = self.___no_cache_key_generator___

        else:
            self._DG_ = mdm_generator
            if cache_key_generator == 'no_cache':
                self._KG_ = self.___no_cache_key_generator___
            else:
                self._KG_ = cache_key_generator

        assert isinstance(correspondence, (list, tuple)), f"put correspondence forms in list please."
        if isinstance(correspondence, tuple): correspondence = list(correspondence)
        self._correspondence_ = correspondence
        self._ndim_ = len(correspondence)

        self._cache_ = dict() # do not use self.___PRIVATE_reset_cache___()
        self.___CT___ = '>CT<'
        self.___NC___ = '>NC<'
        self.___IS_CT___ = False
        self.___IS_NC___ = False
        self.____CT_DG____ = None # the cache for constant data.
        self.___CHECK_repeat_CT___ = True
        self.___CHECK_repeat_CT___ = True
        self.___repeat_CK___ = ''

        self._DO_ = MDM_Do(self)
        self._IS_ = MDM_IS(self)

        self._freeze_self_()

    def ___dict_DG___(self, item):
        return self._dict_mdm_data_[item]

    def ___no_cache_key_generator___(self, item):
        """"""
        assert item in self._iterator_
        return self.___NC___

    @property
    def iterator(self):
        return self._iterator_

    @property
    def correspondence(self):
        return self._correspondence_

    @property
    def ndim(self):
        return self._ndim_

    def __iter__(self):
        """Go through all local basic units."""
        for basic_unit in self._iterator_:
            yield basic_unit

    def __len__(self):
        """How many local basic units."""
        return len(self._iterator_)

    def __contains__(self, item):
        """if item is a local basic unit?"""
        return item in self._iterator_

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

        assert RETURN.ndim == self.ndim, f"ndim of data[{item}] is wrong."
        assert RETURN.__class__.__name__ == 'ndarray', f"I must obtain an ndarray for each basic unit."
        return RETURN

    def __getitem__(self, item):
        """"""
        RETURN = self.___getitem_pre_customizing___(item)
        # customization is after the cache, so we can do whatever customization afterwards.
        # RETURN = self.customize.___PRIVATE_do_execute_customization___(RETURN, item)
        return RETURN

    @property
    def do(self):
        return self._DO_

    @property
    def IS(self):
        return self._IS_

    def __neg__(self):
        """- EWC_SparseMatrix"""
        data_generator = ___NEG___(self)
        RETURN = self.__class__(self.iterator, data_generator, self.correspondence, 'no_cache')
        return RETURN

    def __mul__(self, other):
        """"""
        data_generator = MDM_MUL(self, other)
        RETURN = self.__class__(self.iterator, data_generator, self.correspondence, 'no_cache')
        return RETURN

    def __rmul__(self, other):
        """"""
        data_generator = MDM_MUL(self, other)
        RETURN = self.__class__(self.iterator, data_generator, self.correspondence, 'no_cache')
        return RETURN




if __name__ == '__main__':
    # mpiexec -n 4 python tools/linear_algebra/elementwise_cache/objects/multi_dim_matrix/main.py
    from __init__ import cscg3

    mesh = cscg3.mesh('cuboid', region_layout=[2,2,2])([3,3,3])
    space = cscg3.space('polynomials')((2,2,2))
    FC = cscg3.form(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    velocity = FC('vector', (u,v,w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)

    f1 = FC('1-f', is_hybrid=False)
    u2 = FC('2-f', is_hybrid=False)
    v2 = FC('2-f', is_hybrid=False)

    # f1.TW.func.do.set_func_body_as(velocity)
    # f1.TW.current_time = 0
    # f1.TW.do.push_all_to_instant()
    # f1.discretize()

    MDM = f1.special.cross_product_2f__ip_2f(u2, v2, output='MDM')
    # print(MDM.correspondence)

    # print(MDM.IS.inhomogeneous)
    print(-MDM)