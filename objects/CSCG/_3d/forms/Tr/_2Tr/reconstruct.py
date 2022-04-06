
import numpy as np
from screws.freeze.base import FrozenOnly


class _3dCSCG_2Tr_Reconstruct(FrozenOnly):
    """"""
    def __init__(self, Tr):
        self._Tr_ = Tr
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel=False, i=None):
        """
        Do the reconstruction.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param sigma: A 1d iterable object of floats between -1 and 1.
        :param bool ravel: (`default`:``False``) If we return 1d data?
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :type sigma: list, tuple, numpy.ndarray
        :param i: (`default`:``None``) Do the reconstruction for these
            trace elements. if it is ``None``, then do it for all trace
            elements.
        :type i: int, None, list, tuple
        """
        mesh = self._Tr_.mesh

        # parse indices ------------------------------------------------------------
        if i is None:
            indices = mesh.trace.elements._elements_.keys()
        else:
            if not isinstance(i, (list, tuple)):
                indices = [i,]
            else:
                indices = i

        # ---------------------------------------------------------------------------
        xietasigma, pb = self._Tr_.do.evaluate_basis_at_meshgrid(xi, eta, sigma)
        ii, jj, kk = np.size(xi), np.size(eta), np.size(sigma)
        xyz = dict()
        v = dict()
        for key in indices:
            if key in mesh.trace.elements:
                te = mesh.trace.elements[key]
                side = te.CHARACTERISTIC_side
                ele = te.CHARACTERISTIC_element
                xyz_i = te.coordinate_transformation.mapping(*xietasigma[side],
                                                             from_element=ele,
                                                             side=side)

                g = te.coordinate_transformation.metric(*xietasigma[side])
                prime_cochain = self._Tr_.cochain.local_TEW[key]
                if side in 'NS':
                    vi = np.einsum('i, j, ij -> j', prime_cochain, 1 / np.sqrt(g),
                                   pb[side][0], optimize='greedy')
                elif side in 'WE':
                    vi = np.einsum('i, j, ij -> j', prime_cochain, 1 / np.sqrt(g),
                                   pb[side][0], optimize='greedy')
                elif side in 'BF':
                    vi = np.einsum('i, j, ij -> j', prime_cochain, 1 / np.sqrt(g),
                                   pb[side][0], optimize='greedy')
                else:
                    raise Exception()

                if ravel:
                    xyz[key] = xyz_i
                    v[key] = [vi,]
                else:
                    if side in 'NS':
                        xyz[key] = [xyz_i[m].reshape(jj, kk, order='F') for m in range(3)]
                        v[key] = [vi.reshape((jj, kk), order='F'),]
                    elif side in 'WE':
                        xyz[key] = [xyz_i[m].reshape(ii, kk, order='F') for m in range(3)]
                        v[key] = [vi.reshape((ii, kk), order='F'),]
                    elif side in 'BF':
                        xyz[key] = [xyz_i[m].reshape(ii, jj, order='F') for m in range(3)]
                        v[key] = [vi.reshape((ii, jj), order='F'),]
                    else:
                        raise Exception

        return xyz, v