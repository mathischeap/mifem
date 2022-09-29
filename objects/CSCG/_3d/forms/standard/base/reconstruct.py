# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/29/2022 12:18 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np
from root.config.main import sIze, cOmm, rAnk

___grid_cache_3dCSCG_SF_Reconstruct___ = {
    'grid': None,
    'Xi_Eta_Sigma': None
}


class _3dCSCG_SF_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def discrete_field(self, grid):
        """"""
        if hasattr(self, 'discrete_scalar'):
            return self.discrete_scalar(grid)
        elif hasattr(self, 'discrete_vector'):
            return self.discrete_vector(grid)
        else:
            raise NotImplementedError()

    # noinspection PyTypedDict,PyUnresolvedReferences
    @classmethod
    def ___PRIVATE_distribute_region_wise_meshgrid___(cls, mesh, rst):
        """

        Parameters
        ----------
        rst: dict, list, tuple
            For example:
                rst = [r, s, t], these r, s, t will be used for all regions.
                rst = {'R:R1': [r1, s1, t1], 'R:R2': [r2, s2, t2], ....}, in each
                    region, we use its own r, s, t

            r, s or t can be in [-1,1] or [0,1].

        Returns
        -------

        """
        local_elements_in_regions = mesh.elements.in_regions

        #----------- parse rst -------------------------------------------------------------------1
        if isinstance(rst, (list, tuple)):
            ___ = dict()
            for rn in local_elements_in_regions:
                ___[rn] = rst
            rst = ___

        elif isinstance(rst, dict):
            for rn in rst:
                assert rn in mesh.domain.regions.names, f"rst[{rn}] is not a valid region name."

            for rn in local_elements_in_regions:
                assert rn in rst, f"grid missed for region {rn}."

            ___ = dict()
            for rn in rst:
                if rn in local_elements_in_regions:
                    ___[rn] = rst[rn]
            rst = ___

        else:
            raise Exception(f"rst={rst} wrong.")

        #-------- parse r, s, t -------------------------------------------------------------------1
        for rn in rst:
            r, s, t = rst[rn]
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

            if isinstance(t, (int, float)):
                t = np.array([t, ])
            elif not t.__class__.__name__ == 'ndarray':
                t = np.array(t)
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

            if np.ndim(t) == 1 and min(t) >=0 and max(t) <=1 and np.all(np.diff(t) > 0):
                pass
            elif np.ndim(t) == 1 and min(t) >=-1 and max(t) <=1 and np.all(np.diff(t) > 0):
                t = (t + 1) / 2
            else:
                raise Exception(f"t[{rn}]={t} wrong!")

            rst[rn] = r, s, t

        #------- check cache ----------------------------------------------------------------------1
        if ___grid_cache_3dCSCG_SF_Reconstruct___['grid'] is not None:
            c_rst = ___grid_cache_3dCSCG_SF_Reconstruct___['grid']
            cached = True

            for rn in rst:

                if rn in c_rst:
                    if any([
                        rst[rn][0].shape != c_rst[rn][0].shape,
                        rst[rn][1].shape != c_rst[rn][1].shape,
                        rst[rn][2].shape != c_rst[rn][2].shape
                    ]):
                        cached = False
                        break

                    elif all([
                            np.all(np.abs(rst[rn][0] - c_rst[rn][0]) < 1e-9),
                            np.all(np.abs(rst[rn][1] - c_rst[rn][1]) < 1e-9),
                            np.all(np.abs(rst[rn][2] - c_rst[rn][2]) < 1e-9),
                    ]):
                        pass

                    else:
                        cached = False
                        break

                else:
                    cached = False
                    break

            if cached:
                return ___grid_cache_3dCSCG_SF_Reconstruct___['Xi_Eta_Sigma'], rst
            else:
                pass

        else:
            pass

        #-------- generating RST for region elements ----------------------------------------------1

        Xi_Eta_Sigma_D = dict()
        for rn in rst:
            assert rn in local_elements_in_regions, f"trivial check."
            spacing = mesh.elements.spacing[rn]
            SPr, SPs, SPt = spacing
            L_SPr = len(SPr)
            L_SPs = len(SPs)
            L_SPt = len(SPt)

            RST = list()
            for SP, Len, D in zip(spacing, [L_SPr, L_SPs, L_SPt], rst[rn]):

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

        ___grid_cache_3dCSCG_SF_Reconstruct___['grid']  = rst
        ___grid_cache_3dCSCG_SF_Reconstruct___['Xi_Eta_Sigma'] = Xi_Eta_Sigma_D
        #==========================================================================================1
        return Xi_Eta_Sigma_D, rst

    @classmethod
    def ___PRIVATE_distribute_XYZ_and_VAL___(cls, mesh, xyz, value):
        """"""
        rwPc = mesh.elements.region_wise_prime_core #

        #---------- we distribute xyz, and value to prime cores -------------------------------------1
        XYZ = [dict() for _ in range(sIze)]
        VAL = [dict() for _ in range(sIze)]

        for e in xyz:
            rn = mesh.do.find.region_name_of_element(e)
            prime_core = rwPc[rn]
            XYZ[prime_core][e] = xyz[e]
            VAL[prime_core][e] = value[e]

        XYZ = cOmm.alltoall(XYZ)
        VAL = cOmm.alltoall(VAL)

        #----------- collect xyz and value in the prime cores --------------------------------------1
        _XYZ_ = dict()
        for ___ in XYZ:
            _XYZ_.update(___)
        _VAL_ = dict()
        for ___ in VAL:
            _VAL_.update(___)
        XYZ = _XYZ_
        VAL = _VAL_

        element_global_numbering = dict()
        for rn in rwPc:
            core = rwPc[rn]
            if core == rAnk: # I am the prime core of this region.
                if rn not in element_global_numbering:
                    element_global_numbering[rn] = \
                        mesh.___PRIVATE_generate_element_global_numbering___(number_what=rn)

        return XYZ, VAL, element_global_numbering

    @classmethod
    def ___PRIVATE_prime_region_wise_stack___(cls, mesh, data_dict, dim, rst, element_global_numbering):
        """"""
        _sd_ = dict()

        for Rn in element_global_numbering:
            r, s, t = rst[Rn]
            region_data_shape = len(r), len(s), len(t)

            if dim == 1:
                _3D0 = np.zeros(region_data_shape, dtype='float')

                idk = 0
                for k in range(mesh.elements.layout[Rn][2]):
                    idj = 0
                    for j in range(mesh.elements.layout[Rn][1]):
                        idi = 0
                        for i in range(mesh.elements.layout[Rn][0]):
                            D0 = data_dict[element_global_numbering[Rn][i, j, k]][0]
                            bi, bj, bk = D0.shape
                            _3D0[idi:idi + bi, idj:idj + bj, idk:idk + bk] = D0

                            idi += bi
                        idj += bj
                    idk += bk

                _sd_[Rn] = (_3D0,)

            elif dim == 3:
                _3D0 = np.zeros(region_data_shape, dtype='float')
                _3D1 = np.zeros(region_data_shape, dtype='float')
                _3D2 = np.zeros(region_data_shape, dtype='float')

                idk = 0
                for k in range(mesh.elements.layout[Rn][2]):
                    idj = 0
                    for j in range(mesh.elements.layout[Rn][1]):
                        idi = 0
                        for i in range(mesh.elements.layout[Rn][0]):
                            D0, D1, D2 = data_dict[element_global_numbering[Rn][i, j, k]]
                            bi, bj, bk = D0.shape
                            _3D0[idi:idi+bi, idj:idj+bj, idk:idk+bk] = D0
                            _3D1[idi:idi+bi, idj:idj+bj, idk:idk+bk] = D1
                            _3D2[idi:idi+bi, idj:idj+bj, idk:idk+bk] = D2

                            idi += bi
                        idj += bj
                    idk += bk

                _sd_[Rn] = (_3D0, _3D1, _3D2)

            else:
                raise NotImplementedError()

        return _sd_



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
