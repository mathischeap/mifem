# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/28/2022 12:24 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from screws.freeze.main import FrozenClass
from root.config.main import cOmm, rAnk, mAster_rank


class CSCG_DiscreteField(FrozenClass):
    """Data are region wise. And regions can be distributed in different cores."""

    def __init__(self, mesh, coordinates, values, name, structured=False, grid=None):
        """

        Parameters
        ----------
        mesh
        coordinates
        values
        name
        structured : bool, optional, default: False
            If the data are made from a structured grid?
        grid : {None, list, tuple}, optional
            If the data are made from a structured grid,
            the grid is meshgrid(*grid, indexing='ij').

            So it will only be valid when structured is True.
        """
        self._mesh_ = mesh
        self.standard_properties.name = name

        regions = mesh.domain.regions
        rns = regions.names

        assert isinstance(coordinates, dict), f"put coordinates in a dict pls."
        assert isinstance(values, dict), f"put values in a dict pls."

        _locally_involved_regions_ = set(coordinates.keys())
        assert _locally_involved_regions_ == set(values.keys())
        self._locally_involved_regions_ = tuple(_locally_involved_regions_)

        x_shape = dict()
        for rn in coordinates:
            assert rn in rns, f"coordinates[{rn}] is not a valid region name."
            assert len(coordinates[rn]) == self.ndim, f"coordinates of region {rn} dimension wrong."

            x_shape[rn] = np.shape(coordinates[rn][0])
            yz_ = coordinates[rn][1:]
            for i, _ in enumerate(yz_):
                assert np.shape(_) == x_shape[rn], \
                    f"{i+1} axis shape of coordinates[{rn}] does not match."

        vdim = None
        for rn in values:
            assert rn in rns, f"values[{rn}] is not a valid region name."

            for i, _ in enumerate(values[rn]):
                assert np.shape(_) == x_shape[rn], \
                    f"{i} axis shape {np.shape(_)} of values[{rn}] does not match coordinates shape={x_shape}."

            if vdim is None:
                vdim = len(values[rn])
            else:
                assert vdim == len(values[rn])

        vdim = cOmm.gather(vdim, root=mAster_rank)
        if rAnk == mAster_rank:
            __ck = None
            for _v in vdim:
                if _v is not None:
                    if __ck is None:
                        __ck = _v
                    else:
                        assert __ck == _v
        else:
            __ck = None

        __ck = cOmm.bcast(__ck, root=mAster_rank)

        self._vdim_ = __ck

        self._coordinates_ = coordinates
        self._values_ =values

        if structured: # the data are on a structured grid.
            if grid is not None:
                if not isinstance(grid, dict):
                    _LPD_ = dict()
                    for rn in self.regions:
                        _LPD_[rn] = grid
                    grid = _LPD_

                for rn in self.regions:
                    lsp = grid[rn]

                    for i, _ in enumerate(lsp):
                        assert np.ndim(_) == 1 and np.min(_) >=0 and np.max(_) <= 1, \
                            f"grid[{rn}][{i}] wrong, must be 1d and in [-1,1]."
                        if len(_) > 1:
                            assert np.all(np.diff(_)) > 0, \
                                f"grid[{rn}][{i}] wrong, must be in [-1,1] and increasing."
                        else:
                            pass

                    coo = self.coordinates[rn]
                    val = self.values[rn]

                    SHAPE = tuple([len(_) for _ in lsp])
                    for i, xyz in enumerate(coo):
                        assert xyz.shape == SHAPE, f"{i}th dimension of coordinates dis-match shape={SHAPE}."
                    for i, v in enumerate(val):
                        assert v.shape == SHAPE, f"{i}th dimension of values dis-match shape={SHAPE}."

            else:
                pass
        else:
            pass

        self._structured_ = structured
        self._grid_ = grid

        self._visualize_ = None
        self._portion_ = None
        self._freeze_self_()

    @property
    def mesh(self):
        return self._mesh_

    @property
    def ndim(self):
        return self.mesh.ndim

    @property
    def vdim(self):
        """the dimension of the value. Then vdim == 1, it is a scalar; when vdim == ndim, it
        is a vector; and so on."""
        return self._vdim_
    @property
    def name(self):
        return self.standard_properties.name

    @property
    def coordinates(self):
        return self._coordinates_

    @property
    def values(self):
        return self._values_

    @property
    def structured(self):
        """Return True is the data are on a structured grid else return False."""
        return self._structured_

    @property
    def grid(self):
        """If `self.structured`, the grid is `self.grid`"""
        return self._grid_

    @property
    def regions(self):
        """The locally involved regions. So data are actually stored in different regions distributed
        in different cores.

        The locally involved regions will be according to the region-wise-prime-core property of
        the mesh. So it is fixed for a given mesh.
        """
        return self._locally_involved_regions_

    @property
    def visualize(self):
        """Try to visualize the data."""
        return self._visualize_

    @property
    def portion(self):
        """A portion of the data."""
        return self._portion_


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
