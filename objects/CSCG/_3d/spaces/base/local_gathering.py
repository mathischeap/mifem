# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/2/2022 5:31 PM
"""
from components.freeze.main import FrozenOnly
import numpy as np


class _3dCSCG_Space_LocalGathering(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()

    @property
    def _3dCSCG_0Trace(self):
        px, py, pz = self._space_.p
        fully_local_numbering_3d = - np.ones((px+1, py+1, pz+1), dtype=int)
        cn = 0
        fully_local_numbering_3d[0, :, :] = np.arange(cn, cn+(py+1)*(pz+1)).reshape((py+1, pz+1), order='F')
        cn += (py + 1) * (pz + 1)
        fully_local_numbering_3d[-1, :, :] = np.arange(cn, cn+(py+1)*(pz+1)).reshape((py+1, pz+1), order='F')
        cn += (py + 1) * (pz + 1)
        fully_local_numbering_3d[1:-1, 0, :] = np.arange(cn, cn+(px-1)*(pz+1)).reshape((px-1, pz+1), order='F')
        cn += (px - 1) * (pz + 1)
        fully_local_numbering_3d[1:-1, -1, :] = np.arange(cn, cn+(px-1)*(pz+1)).reshape((px-1, pz+1), order='F')
        cn += (px - 1) * (pz + 1)
        fully_local_numbering_3d[1:-1, 1:-1, 0] = np.arange(cn, cn+(px-1)*(py-1)).reshape((px-1, py-1), order='F')
        cn += (px - 1) * (py - 1)
        fully_local_numbering_3d[1:-1, 1:-1, -1] = np.arange(cn, cn+(px-1)*(py-1)).reshape((px-1, py-1), order='F')

        local_gathering = {
            'N': fully_local_numbering_3d[0, :, :].ravel('F'),
            'S': fully_local_numbering_3d[-1, :, :].ravel('F'),
            'W': fully_local_numbering_3d[:, 0, :].ravel('F'),
            'E': fully_local_numbering_3d[:, -1, :].ravel('F'),
            'B': fully_local_numbering_3d[:, :, 0].ravel('F'),
            'F': fully_local_numbering_3d[:, :, -1].ravel('F'),
        }

        for side in local_gathering:
            assert -1 not in local_gathering[side], f"local numbering for side {side} not right."

        return local_gathering

    @property
    def _3dCSCG_1Trace(self):
        px, py, pz = self._space_.p
        dx_3d = - np.ones((px, py+1, pz+1), dtype=int)
        dy_3d = - np.ones((px+1, py, pz+1), dtype=int)
        dz_3d = - np.ones((px+1, py+1, pz), dtype=int)

        cn = 0

        dx_3d[:, 0, :] = np.arange(cn, cn + px * (pz + 1)).reshape((px, pz + 1), order='F')
        cn += px * (pz + 1)
        dx_3d[:, -1, :] = np.arange(cn, cn + px * (pz + 1)).reshape((px, pz + 1), order='F')
        cn += px * (pz + 1)
        dx_3d[:, 1:-1, 0] = np.arange(cn, cn + px * (py - 1)).reshape((px, py - 1), order='F')
        cn += px * (py - 1)
        dx_3d[:, 1:-1, -1] = np.arange(cn, cn + px * (py - 1)).reshape((px, py - 1), order='F')
        cn += px * (py - 1)

        dy_3d[0, :, :] = np.arange(cn, cn + py * (pz + 1)).reshape((py, pz + 1), order='F')
        cn += py * (pz + 1)
        dy_3d[-1, :, :] = np.arange(cn, cn + py * (pz + 1)).reshape((py, pz + 1), order='F')
        cn += py * (pz + 1)
        dy_3d[1:-1, :, 0] = np.arange(cn, cn + py * (px - 1)).reshape((px-1, py), order='F')
        cn += py * (px - 1)
        dy_3d[1:-1, :, -1] = np.arange(cn, cn + py * (px - 1)).reshape((px-1, py), order='F')
        cn += py * (px - 1)

        dz_3d[0, :, :] = np.arange(cn, cn + pz * (py + 1)).reshape((py + 1, pz), order='F')
        cn += pz * (py + 1)
        dz_3d[-1, :, :] = np.arange(cn, cn + pz * (py + 1)).reshape((py + 1, pz), order='F')
        cn += pz * (py + 1)
        dz_3d[1:-1, 0, :] = np.arange(cn, cn + pz * (px - 1)).reshape((px - 1, pz), order='F')
        cn += pz * (px - 1)
        dz_3d[1:-1, -1, :] = np.arange(cn, cn + pz * (px - 1)).reshape((px - 1, pz), order='F')

        local_gathering = {
            'N': np.concatenate([dy_3d[0, :, :].ravel('F'), dz_3d[0, :, :].ravel('F')]),
            'S': np.concatenate([dy_3d[-1, :, :].ravel('F'), dz_3d[-1, :, :].ravel('F')]),
            'W': np.concatenate([dx_3d[:, 0, :].ravel('F'), dz_3d[:, 0, :].ravel('F')]),
            'E': np.concatenate([dx_3d[:, -1, :].ravel('F'), dz_3d[:, -1, :].ravel('F')]),
            'B': np.concatenate([dx_3d[:, :, 0].ravel('F'), dy_3d[:, :, 0].ravel('F')]),
            'F': np.concatenate([dx_3d[:, :, -1].ravel('F'), dy_3d[:, :, -1].ravel('F')]),
        }

        for side in local_gathering:
            assert -1 not in local_gathering[side], f"local numbering for side {side} not right."

        return local_gathering

    @property
    def _3dCSCG_2Trace(self):
        px, py, pz = self._space_.p
        local_gathering = dict()
        cn = 0
        local_gathering['N'] = [i for i in range(cn, cn + py * pz)]
        cn += py * pz
        local_gathering['S'] = [i for i in range(cn, cn + py * pz)]
        cn += py * pz
        local_gathering['W'] = [i for i in range(cn, cn + px * pz)]
        cn += px * pz
        local_gathering['E'] = [i for i in range(cn, cn + px * pz)]
        cn += px * pz
        local_gathering['B'] = [i for i in range(cn, cn + px * py)]
        cn += px * py
        local_gathering['F'] = [i for i in range(cn, cn + px * py)]

        for side in local_gathering:
            assert -1 not in local_gathering[side], f"local numbering for side {side} not right."

        return local_gathering

    @property
    def _3dCSCG_0LocalTrace(self):
        return self._3dCSCG_0Trace

    @property
    def _3dCSCG_1LocalTrace(self):
        return self._3dCSCG_1Trace

    @property
    def _3dCSCG_2LocalTrace(self):
        return self._3dCSCG_2Trace
