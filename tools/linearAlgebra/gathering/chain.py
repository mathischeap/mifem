# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/09 10:04 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix
from tools.linearAlgebra.gathering.irregular.ir_chain_matrix.main import iR_Chain_Gathering_Matrix


___GLOBAL_CGM_CACHE___ = {'key': list(), 'cache': list()} # to cache global CGM.


class GatheringMatrixChaining(FrozenOnly):
    """"""

    def __init__(self, *forms):
        """"""
        GM_names = list()
        GMs = list()
        for f in forms:
            if f.__class__.__name__ in ('Gathering_Matrix', 'iR_Gathering_Matrix'):
                GM = f
            else:
                GM = f.numbering.gathering
            GMs.append(GM)
            GM_names.append(GM.__class__.__name__)

        self._GMs_ = GMs
        self._GM_names_ = GM_names

        self._freeze_self_()


    def __call__(self, chain_method='silly'):
        """"""
        cache_key = ''.join([_.stamp for _ in self._GMs_]) + chain_method
        if cache_key in ___GLOBAL_CGM_CACHE___['key']:
            index = ___GLOBAL_CGM_CACHE___['key'].index(cache_key)
            return ___GLOBAL_CGM_CACHE___['cache'][index]

        else:
            if all([_ == 'Gathering_Matrix' for _ in self._GM_names_]):

                CGM = Chain_Gathering_Matrix(self._GMs_, chain_method=chain_method)

            elif any([_ == 'iR_Gathering_Matrix' for _ in self._GM_names_]):

                CGM = iR_Chain_Gathering_Matrix(self._GMs_, chain_method=chain_method)

            else:
                raise NotImplementedError()

            ___GLOBAL_CGM_CACHE___['key'].append(cache_key)
            ___GLOBAL_CGM_CACHE___['cache'].append(CGM)

            if len(___GLOBAL_CGM_CACHE___['key']) > 3: # only cache 3 CGM
                ___GLOBAL_CGM_CACHE___['key'] = ___GLOBAL_CGM_CACHE___['key'][-3:]
                ___GLOBAL_CGM_CACHE___['cache'] = ___GLOBAL_CGM_CACHE___['cache'][-3:]
            else:
                pass

            return CGM


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
