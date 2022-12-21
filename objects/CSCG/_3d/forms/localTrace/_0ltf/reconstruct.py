# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:55 PM
"""
from components.freeze.main import FrozenOnly
import numpy as np


class _3dCSCG_0LocalTrace_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel=False, element_range=None):
        """

        Parameters
        ----------
        xi
        eta
        sigma
        ravel
        element_range : mesh element number

        Returns
        -------

        """
        mesh = self._ltf_.mesh

        if element_range is None:
            indices = mesh.elements._elements_.keys()
        elif isinstance(element_range, int):
            indices = [element_range, ]
        else:
            indices = element_range

        xietasigma, pb = self._ltf_.do.evaluate_basis_at_meshgrid(xi, eta, sigma)
        ii, jj, kk = np.size(xi), np.size(eta), np.size(sigma)
        xyz = dict()
        v = dict()

        Tmap = mesh.trace.elements.map

        for i in indices:
            if i in mesh.elements:

                tes = Tmap[i]
                xyz_i = list()
                v_i = list()
                for te_number, side in zip(tes, 'NSWEBF'):
                    te = mesh.trace.elements[te_number]

                    if i in self._ltf_.cochain.local:
                        side_prime_cochain = \
                            self._ltf_.cochain.___components_of_cochain_on_element_side___(
                                i, side
                            )
                    else:
                        if i in self._ltf_.cochain.local_ESW:
                            esw_i = self._ltf_.cochain.local_ESW[i]
                            if side in esw_i:
                                side_prime_cochain = esw_i[side]
                            else:
                                side_prime_cochain = None
                        else:
                            side_prime_cochain = None

                    if side_prime_cochain is not None:
                        xyz_i.append(
                            te.coordinate_transformation.mapping(
                                *xietasigma[side], from_element=i, side=side
                            )
                        )
                        v_i.append(
                            np.einsum(
                                'i, ij -> j',
                                side_prime_cochain,
                                pb[side][0],
                                optimize='greedy'
                            )
                        )

                    else:
                        xyz_i.append(None)
                        v_i.append(None)

                if all([_ is None for _ in v_i]):
                    continue
                else:
                    pass

                XYZ = dict()
                V = dict()

                if ravel:
                    for j, side in enumerate('NSWEBF'):
                        XYZ[side] = xyz_i[j]
                        V[side] = [v_i[j], ]

                else:
                    for j, side in enumerate('NSWEBF'):
                        xyz_ = xyz_i[j]
                        v_ = v_i[j]

                        if xyz_ is None:
                            pass
                        else:
                            if side in 'NS':
                                xyz_ = [xyz_[_].reshape((jj, kk), order='F') for _ in range(3)]
                                v_ = [v_.reshape((jj, kk), order='F'), ]
                            elif side in 'WE':
                                xyz_ = [xyz_[_].reshape((ii, kk), order='F') for _ in range(3)]
                                v_ = [v_.reshape((ii, kk), order='F'), ]
                            elif side in 'BF':
                                xyz_ = [xyz_[_].reshape((ii, jj), order='F') for _ in range(3)]
                                v_ = [v_.reshape((ii, jj), order='F'), ]
                            else:
                                raise Exception

                        XYZ[side] = xyz_
                        V[side] = v_

                xyz[i] = XYZ
                v[i] = V

        return xyz, v

    def ___PrLT_make_reconstruction_matrix_on_grid___(self, xi, et, sg, element_range=None):
        """

        Parameters
        ----------
        xi
        et
        sg
        element_range

        Returns
        -------
        RD : dict
            A dict contains mesh-element-wise reconstruction matrix.
        RD : dict
            A dict (keys are mesh-element-wise) of dictionaries which
            are element-side-wise reconstruction matrix.

        """
        mesh = self._ltf_.mesh

        if element_range is None:
            indices = mesh.elements._elements_.keys()

        elif isinstance(element_range, int):
            indices = [element_range, ]

        elif isinstance(element_range, list):  # a list of positions.

            # a list of positions, we assume they are element_sides or something else!

            indices = dict()
            for pos in element_range:
                if isinstance(pos, str) and len(pos) > 1 and pos[0].isnumeric() and not pos[-1].isnumeric():
                    element = int(pos[:-1])
                    side = pos[-1]

                    if element not in indices:
                        indices[element] = list()
                    else:
                        pass

                    if side in indices[element]:
                        pass
                    else:
                        indices[element].append(side)

                else:
                    raise NotImplementedError(f"position {pos} is not not understandable.")

        else:
            indices = element_range

        xietasigma, pb = self._ltf_.do.evaluate_basis_at_meshgrid(xi, et, sg)
        RD_sides = dict()

        M = list()
        for side in 'NSWEBF':
            M.append(pb[side][0].T)

        rdi_sides = {
            'N': M[0],
            'S': M[1],
            'W': M[2],
            'E': M[3],
            'B': M[4],
            'F': M[5],
        }

        if not isinstance(indices, dict):
            for i in indices:
                if i in mesh.elements:
                    RD_sides[i] = rdi_sides

            return RD_sides

        else:  # we get a dictionary indices whose keys are then mesh-elements and values are mesh-element-positions.

            for i in indices:
                if i in mesh.elements:
                    RD_sides[i] = dict()
                    for side in indices[i]:
                        RD_sides[i][side] = rdi_sides[side]
                else:
                    pass

            return RD_sides
