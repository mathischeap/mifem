# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/30 2:38 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np
from root.config.main import sIze, cOmm, rAnk

___grid_cache_2dCSCG_SF_Reconstruct___ = {
    'grid': None,
    'Xi_Eta_Sigma': None
}


class _2dCSCG_SF_ReconstructBase(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def discrete_field(self, grid):
        """"""
        if hasattr(self, 'discrete_scalar'):
            return self.discrete_scalar(grid)
        elif hasattr(self, 'discrete_vector'):
            return self.discrete_vector(grid)
        else:
            raise NotImplementedError()

    # noinspection PyTypedDict
    @classmethod
    def ___PRIVATE_distribute_region_wise_meshgrid___(cls, mesh, rs):
        """

        Parameters
        ----------
        rs: dict, list, tuple
            For example:
                rst = [r, s], these r, s will be used for all regions.
                rst = {'R:R1': [r1, s1], 'R:R2': [r2, s2], ....}, in each
                    region, we use its own r, s.

            r or s can be in [-1,1] or [0,1].

        Returns
        -------

        """

        local_elements_in_regions = mesh.elements.in_regions

        #----------- parse rst -------------------------------------------------------------------1
        if isinstance(rs, (list, tuple)):
            ___ = dict()
            for rn in local_elements_in_regions:
                # noinspection PyUnresolvedReferences
                ___[rn] = rs
            rs = ___

        elif isinstance(rs, dict):
            for rn in rs:
                assert rn in mesh.domain.regions.names, f"rs[{rn}] is not a valid region name."

            for rn in local_elements_in_regions:
                assert rn in rs, f"grid missed for region {rn}."

            ___ = dict()
            for rn in rs:
                if rn in local_elements_in_regions:
                    ___[rn] = rs[rn]
            rs = ___

        else:
            raise Exception(f"grid={rs} wrong.")

        #-------- parse r, s further --------------------------------------------------------------1
        for rn in rs:
            r, s = rs[rn]
            if isinstance(r, (int, float)):
                r = np.array([r, ])
            elif not r.__class__.__name__ == 'ndarray':
                r = np.array(r)
            else:
                pass

            if isinstance(s, (int, float)):
                s = np.array([s, ])
            elif not s.__class__.__name__ == 'ndarray':
                s = np.array(s)
            else:
                pass

            if np.ndim(r) == 1 and min(r) >=0 and max(r) <=1 and np.all(np.diff(r) > 0):
                pass
            elif np.ndim(r) == 1 and min(r) >=-1 and max(r) <=1 and np.all(np.diff(r) > 0):
                r = (r + 1) / 2
            else:
                raise Exception(f"r[{rn}]={r} wrong!")

            if np.ndim(s) == 1 and min(s) >=0 and max(s) <=1 and np.all(np.diff(s) > 0):
                pass
            elif np.ndim(s) == 1 and min(s) >=-1 and max(s) <=1 and np.all(np.diff(s) > 0):
                s = (s + 1) / 2
            else:
                raise Exception(f"s[{rn}]={s} wrong!")

            rs[rn] = r, s

        #------- check cache ----------------------------------------------------------------------1
        if ___grid_cache_2dCSCG_SF_Reconstruct___['grid'] is not None:
            c_rs = ___grid_cache_2dCSCG_SF_Reconstruct___['grid']
            cached = True
            for rn in rs:
                if rn in c_rs:
                    # noinspection PyUnresolvedReferences
                    if any([
                        rs[rn][0].shape != c_rs[rn][0].shape,
                        rs[rn][1].shape != c_rs[rn][1].shape,
                    ]):
                        cached = False
                        break
                    elif all([
                            np.all(np.abs(rs[rn][0] - c_rs[rn][0]) < 1e-9),
                            np.all(np.abs(rs[rn][1] - c_rs[rn][1]) < 1e-9),
                    ]):
                        pass
                    else:
                        cached = False
                        break

                else:
                    cached = False
                    break
            if cached:
                return ___grid_cache_2dCSCG_SF_Reconstruct___['Xi_Eta_Sigma'], rs
            else:
                pass
        else:
            pass

        #-------- generating RST for region elements ----------------------------------------------1

        Xi_Eta_Sigma_D = dict()
        for rn in rs:
            assert rn in local_elements_in_regions, f"trivial check."
            spacing = mesh.elements.spacing[rn]
            SPr, SPs = spacing
            L_SPr = len(SPr)
            L_SPs = len(SPs)

            RST = list()
            for SP, Len, D in zip(spacing, [L_SPr, L_SPs], rs[rn]):

                Segments = [list() for _ in range(Len-1)]

                base_i = 0
                for _ in D:
                    for i, low_bound in enumerate(SP[base_i:-1]):
                        I = base_i + i
                        up_bound = SP[I +1]
                        if low_bound <= _ < up_bound:
                            Segments[I].append(_)
                            base_i = I
                            break

                if D[-1] == 1:
                    if Segments[-1] == list():
                        Segments[-1].append(1)
                    else:
                        if Segments[-1][-1] != 1.0:
                            Segments[-1].append(1)
                        else:
                            pass

                for i, seg in enumerate(Segments):
                    Seg_Len = SP[i+1] - SP[i]
                    ARR = (np.array(seg) - SP[i]) * 2 / Seg_Len - 1
                    Segments[i] = ARR

                RST.append(Segments)

            Xi_Eta_Sigma_D[rn] = RST

        ___grid_cache_2dCSCG_SF_Reconstruct___['grid']  = rs
        ___grid_cache_2dCSCG_SF_Reconstruct___['Xi_Eta_Sigma'] = Xi_Eta_Sigma_D
        #======================================================================================
        return Xi_Eta_Sigma_D, rs

    @classmethod
    def ___PRIVATE_distribute_XY_and_VAL___(cls, mesh, xy, value):
        """"""
        rwPc = mesh.elements.region_wise_prime_core #

        #---------- we distribute xyz, and value to prime cores -------------------------------------1
        XY = [dict() for _ in range(sIze)]
        VAL = [dict() for _ in range(sIze)]

        for e in xy:
            rn = mesh.do.find.region_name_of_element(e)
            prime_core = rwPc[rn]
            XY[prime_core][e] = xy[e]
            VAL[prime_core][e] = value[e]

        XY = cOmm.alltoall(XY)
        VAL = cOmm.alltoall(VAL)

        #----------- collect xyz and value in the prime cores --------------------------------------1
        _XY_ = dict()
        for ___ in XY:
            _XY_.update(___)
        _VAL_ = dict()
        for ___ in VAL:
            _VAL_.update(___)
        XY = _XY_
        VAL = _VAL_

        element_global_numbering = dict()
        for rn in rwPc:
            core = rwPc[rn]
            if core == rAnk: # I am the prime core of this region.
                if rn not in element_global_numbering:
                    element_global_numbering[rn] = \
                        mesh.___PRIVATE_generate_element_global_numbering___(number_what=rn)

        return XY, VAL, element_global_numbering

    @classmethod
    def ___PRIVATE_prime_region_wise_stack___(cls, mesh, data_dict, dim, rs, element_global_numbering):
        """"""
        _sd_ = dict()

        for Rn in element_global_numbering:
            r, s = rs[Rn]
            region_data_shape = len(r), len(s)

            if dim == 1:
                _3D0 = np.zeros(region_data_shape, dtype='float')

                idj = 0
                for j in range(mesh.elements.layout[Rn][1]):
                    idi = 0
                    for i in range(mesh.elements.layout[Rn][0]):
                        D0 = data_dict[element_global_numbering[Rn][i, j]][0]
                        bi, bj = D0.shape
                        _3D0[idi:idi + bi, idj:idj + bj] = D0

                        idi += bi
                    idj += bj

                _sd_[Rn] = (_3D0,)

            elif dim == 2:
                _3D0 = np.zeros(region_data_shape, dtype='float')
                _3D1 = np.zeros(region_data_shape, dtype='float')

                idj = 0
                for j in range(mesh.elements.layout[Rn][1]):
                    idi = 0
                    for i in range(mesh.elements.layout[Rn][0]):
                        D0, D1 = data_dict[element_global_numbering[Rn][i, j]]
                        bi, bj = D0.shape
                        _3D0[idi:idi+bi, idj:idj+bj] = D0
                        _3D1[idi:idi+bi, idj:idj+bj] = D1

                        idi += bi
                    idj += bj

                _sd_[Rn] = (_3D0, _3D1)

            else:
                raise NotImplementedError()

        return _sd_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
