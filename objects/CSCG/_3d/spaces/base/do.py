from screws.freeze.main import FrozenOnly

import numpy as np
from screws.exceptions import DimensionError
from screws.decorators.accepts import memoize1

class _3dCSCG_space_do(FrozenOnly):
    """"""
    def __init__(self, space):
        self._space_ = space
        self.ndim = space.ndim
        self.basises = space.basises
        self._freeze_self_()

    def evaluate_edge_basis_at_meshgrid(self, k, *domain):
        """ """
        assert 0 <= k <= self.ndim-1, " <Polynomials> : k={} is wrong.".format(k)
        assert len(domain) == self.ndim, \
            " <Polynomials> : domain shape={} wrong.".format(np.shape(domain))
        for i in range(self.ndim):
            assert domain[i].__class__.__name__ in ('list', 'ndarray'), \
                " <Polynomials> : domain[{}].type={} is wrong.".format(i, domain[i].__class__.__name__)
            assert np.ndim(domain[i]) == 1, \
                " <Polynomials> : ndim(domain[{}])={} is wrong.".format(i, np.ndim(domain[i]))
            if np.size(domain[i]) > 1:
                assert np.all(np.diff(domain[i]) > 0) and np.max(domain[i]) <= 1 and np.min(domain[i]) >= -1, \
                    " <Polynomials> : domain[i]={} wrong, need to be increasing and bounded in [-1, 1].".format(
                            domain[i])
            else:
                pass

        if self.ndim == 3:
            if k == 0:
                bf_NS = self.basises[0].node_basis(x=domain[0])
                bf_WE = self.basises[1].node_basis(x=domain[1])
                bf_BF = self.basises[2].node_basis(x=domain[2])
            elif k == 1:
                bf_NS = self.basises[0].edge_basis(x=domain[0])
                bf_WE = self.basises[1].edge_basis(x=domain[1])
                bf_BF = self.basises[2].edge_basis(x=domain[2])
            else:
                raise Exception()
        else:
            raise Exception()

        bNS = (bf_NS,)
        bWE = (bf_WE,)
        bBF = (bf_BF,)

        _basis_ = {'WB': bNS,
                   'EB': bNS,
                   'WF': bNS,
                   'EF': bNS,
                   'NB': bWE,
                   'SB': bWE,
                   'NF': bWE,
                   'SF': bWE,
                   'NW': bBF,
                   'SW': bBF,
                   'NE': bBF,
                   'SE': bBF}

        return _basis_

    def evaluate_form_basis_at_meshgrid(self, k, *domain, compute_xietasigma=True):
        """
        Parameters
        ---------
        k : int
        domain : tuple
            The domain we are going to evaluate basis. Notice that here we only
            accept tuple. So even for 1-D Polynomials, we have to put xi in a
            tuple: (xi,).
        compute_xietasigma :
        """
        assert 0 <= k <= self.ndim, " <Polynomials> : k={} is wrong.".format(k)
        assert len(domain) == self.ndim, \
            " <Polynomials> : domain shape={} wrong.".format(np.shape(domain))
        for i in range(self.ndim):
            assert domain[i].__class__.__name__ in ('list', 'ndarray'), \
                " <Polynomials> : domain[{}].type={} is wrong.".format(
                    i, domain[i].__class__.__name__)
            assert np.ndim(domain[i]) == 1, \
                " <Polynomials> : ndim(domain[{}])={} is wrong.".format(
                    i, np.ndim(domain[i]))
            if np.size(domain[i]) > 1:
                assert np.all(np.diff(domain[i]) > 0) and np.max(domain[i]) <= 1 and np.min(
                    domain[i]) >= -1, \
                    " <Polynomials> : domain[i]={} wrong, need to be increasing and bounded in [-1, 1].".format(
                        domain[i])
            else:
                pass

        if self.ndim == 3:
            if compute_xietasigma:
                _xietasigma_ = self.___PRIVATE_evaluate_basis_functions_meshgrid_3d___(domain)
            if k == 0:
                bf_xi = self.basises[0].node_basis(x=domain[0])
                bf_et = self.basises[1].node_basis(x=domain[1])
                bf_si = self.basises[2].node_basis(x=domain[2])
                bf = np.kron(np.kron(bf_si, bf_et), bf_xi)

                _basis_ = (bf,)

            elif k == 1:
                ed_xi = self.basises[0].edge_basis(x=domain[0])
                lb_et = self.basises[1].node_basis(x=domain[1])
                lb_si = self.basises[2].node_basis(x=domain[2])
                bf_edge_dxi = np.kron(np.kron(lb_si, lb_et), ed_xi)

                lb_xi = self.basises[0].node_basis(x=domain[0])
                ed_et = self.basises[1].edge_basis(x=domain[1])
                lb_si = self.basises[2].node_basis(x=domain[2])
                bf_edge_det = np.kron(np.kron(lb_si, ed_et), lb_xi)

                lb_xi = self.basises[0].node_basis(x=domain[0])
                lb_et = self.basises[1].node_basis(x=domain[1])
                ed_si = self.basises[2].edge_basis(x=domain[2])
                bf_edge_dsi = np.kron(np.kron(ed_si, lb_et), lb_xi)

                _basis_ = (bf_edge_dxi, bf_edge_det, bf_edge_dsi)

            elif k == 2:
                lb_xi = self.basises[0].node_basis(x=domain[0])
                ed_et = self.basises[1].edge_basis(x=domain[1])
                ed_si = self.basises[2].edge_basis(x=domain[2])
                bf_face_det_dsi = np.kron(np.kron(ed_si, ed_et), lb_xi)

                ed_xi = self.basises[0].edge_basis(x=domain[0])
                lb_et = self.basises[1].node_basis(x=domain[1])
                ed_si = self.basises[2].edge_basis(x=domain[2])
                bf_face_dsi_dxi = np.kron(np.kron(ed_si, lb_et), ed_xi)

                ed_xi = self.basises[0].edge_basis(x=domain[0])
                ed_et = self.basises[1].edge_basis(x=domain[1])
                lb_si = self.basises[2].node_basis(x=domain[2])
                bf_face_dxi_det = np.kron(np.kron(lb_si, ed_et), ed_xi)

                _basis_ = (bf_face_det_dsi, bf_face_dsi_dxi, bf_face_dxi_det)

            elif k == 3:
                bf_xi = self.basises[0].edge_basis(x=domain[0])
                bf_et = self.basises[1].edge_basis(x=domain[1])
                bf_si = self.basises[2].edge_basis(x=domain[2])
                bf = np.kron(np.kron(bf_si, bf_et), bf_xi)

                _basis_ = (bf,)
            else:
                raise DimensionError()
        else:
            raise DimensionError()

        if compute_xietasigma:
            # noinspection PyUnboundLocalVariable
            return _xietasigma_, _basis_
        else:
            return None, _basis_

    @staticmethod
    def ___PRIVATE_evaluate_basis_functions_meshgrid_3d___(domain):
        xi, eta, sigma = np.meshgrid(*domain, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        sigma = sigma.ravel('F')
        return xi, eta, sigma

    @memoize1
    def evaluate_form_basis_at_quadrature(
            self, k, quad_degree, quad_type=None, compute_xietasigma=True):
        """"""
        quad_nodes, quad_weights, quad_weights_ravel = \
            self._space_.___PRIVATE_do_evaluate_quadrature___(quad_degree, quad_type=quad_type)
        _xietasigma_, _basis_ = \
            self.evaluate_form_basis_at_meshgrid(
                k, *quad_nodes, compute_xietasigma=compute_xietasigma)
        return _xietasigma_, _basis_, quad_weights, quad_weights_ravel

    def evaluate_trace_basis_at_meshgrid(self, k, *domain, compute_xietasigma=True):
        """ """
        assert 0 <= k <= self.ndim-1, " <Polynomials> : k={} is wrong.".format(k)
        assert len(domain) == self.ndim, \
            " <Polynomials> : domain shape={} wrong.".format(np.shape(domain))
        for i in range(self.ndim):
            assert domain[i].__class__.__name__ in ('list', 'ndarray'), \
                " <Polynomials> : domain[{}].type={} is wrong.".format(i, domain[i].__class__.__name__)
            assert np.ndim(domain[i]) == 1, \
                " <Polynomials> : ndim(domain[{}])={} is wrong.".format(i, np.ndim(domain[i]))
            if np.size(domain[i]) > 1:
                assert np.all(np.diff(domain[i]) > 0) and np.max(domain[i]) <= 1 and np.min(domain[i]) >= -1, \
                    " <Polynomials> : domain[i]={} wrong, need to be increasing and bounded in [-1, 1].".format(
                            domain[i])
            else:
                pass
        if self.ndim == 3:

            if compute_xietasigma:
                # noinspection PyUnresolvedReferences
                _meshgrid_ravel_ = \
                    {'N': [item.ravel('F') for item in np.meshgrid([-1], domain[1], domain[2], indexing='ij')],
                     'S': [item.ravel('F') for item in np.meshgrid([+1], domain[1], domain[2], indexing='ij')],
                     'W': [item.ravel('F') for item in np.meshgrid(domain[0], [-1], domain[2], indexing='ij')],
                     'E': [item.ravel('F') for item in np.meshgrid(domain[0], [+1], domain[2], indexing='ij')],
                     'B': [item.ravel('F') for item in np.meshgrid(domain[0], domain[1], [-1], indexing='ij')],
                     'F': [item.ravel('F') for item in np.meshgrid(domain[0], domain[1], [+1], indexing='ij')]
                     }
            else:
                _meshgrid_ravel_ = None

            _basis_ = dict()
            if k == 0:
                bf_xi = self.basises[0].node_basis(x=domain[0])
                bf_et = self.basises[1].node_basis(x=domain[1])
                bf_si = self.basises[2].node_basis(x=domain[2])
                _basis_['N'] = (np.kron(bf_si, bf_et),)
                _basis_['W'] = (np.kron(bf_si, bf_xi),)
                _basis_['B'] = (np.kron(bf_et, bf_xi),)
            elif k == 1:
                bfn_xi = self.basises[0].node_basis(x=domain[0])
                bfn_et = self.basises[1].node_basis(x=domain[1])
                bfn_si = self.basises[2].node_basis(x=domain[2])
                bfe_xi = self.basises[0].edge_basis(x=domain[0])
                bfe_et = self.basises[1].edge_basis(x=domain[1])
                bfe_si = self.basises[2].edge_basis(x=domain[2])
                _basis_['N'] = (np.kron(bfn_si, bfe_et), np.kron(bfe_si, bfn_et))
                _basis_['W'] = (np.kron(bfn_si, bfe_xi), np.kron(bfe_si, bfn_xi))
                _basis_['B'] = (np.kron(bfn_et, bfe_xi), np.kron(bfe_et, bfn_xi))
            elif k == 2:
                bf_xi = self.basises[0].edge_basis(x=domain[0])
                bf_et = self.basises[1].edge_basis(x=domain[1])
                bf_si = self.basises[2].edge_basis(x=domain[2])
                _basis_['N'] = (np.kron(bf_si, bf_et),)
                _basis_['W'] = (np.kron(bf_si, bf_xi),)
                _basis_['B'] = (np.kron(bf_et, bf_xi),)
            else:
                raise DimensionError()
        else:
            raise DimensionError(
                " <Polynomials> : evaluate_basis_functions for ndim={} not coded".format(self.ndim))
        _basis_['S'] = _basis_['N']
        _basis_['E'] = _basis_['W']
        _basis_['F'] = _basis_['B']
        return _meshgrid_ravel_, _basis_