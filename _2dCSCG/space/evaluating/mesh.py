
import numpy as np
from SCREWS.exceptions import DimensionError

# noinspection PyUnresolvedReferences
class EvaluatingMeshBasis:
    def DO_evaluate_form_basis_at_meshgrid(self, k, *domain, orientation=None, compute_xietasigma=True):
        """
        Parameters
        ---------
        k : int
        domain : tuple
            The domain we are going to evaluate basis. Notice that here we only
            accept tuple. So even for 1-D Polynomials, we have to put xi in a
            tuple: (xi,).
        orientation : In 2d CSCG mesh, orientation really matters!
        compute_xietasigma :
        """
        assert 0 <= k <= self.ndim, " <2dCSCG Space> : k={} is wrong.".format(k)
        assert len(domain) == self.ndim == 2, \
            " <2dCSCG Space> : domain shape={} wrong.".format(np.shape(domain))
        for i in range(self.ndim):
            assert domain[i].__class__.__name__ in ('list', 'ndarray'), \
                " <2dCSCG Space> : domain[{}].type={} is wrong.".format(i, domain[i].__class__.__name__)
            assert np.ndim(domain[i]) == 1, \
                " <2dCSCG Space> : ndim(domain[{}])={} is wrong.".format(i, np.ndim(domain[i]))
            if np.size(domain[i]) > 1:
                assert np.all(np.diff(domain[i]) > 0) and \
                       np.max(domain[i]) <= 1 and \
                       np.min(domain[i]) >= -1, \
                    " <2dCSCG Space> : domain[i]={} wrong, need to be " \
                    "increasing and bounded in [-1, 1].".format( domain[i])
            else:
                pass

        # ... 0 ...
        if k == 0:
            bf_xi = self.basises[0].node_basis(x=domain[0])
            bf_et = self.basises[1].node_basis(x=domain[1])
            bf = np.kron(bf_et, bf_xi)
            _basis_ = (bf,)
        # ... 1 ...
        elif k == 1:
            if orientation == 'inner':
                ed_xi = self.basises[0].edge_basis(x=domain[0])
                lb_et = self.basises[1].node_basis(x=domain[1])
                bf_edge_dxi = np.kron(lb_et, ed_xi)
                lb_xi = self.basises[0].node_basis(x=domain[0])
                ed_et = self.basises[1].edge_basis(x=domain[1])
                bf_edge_det = np.kron(ed_et, lb_xi)
                _basis_ = (bf_edge_dxi, bf_edge_det)
            elif orientation == 'outer':
                lb_xi = self.basises[0].node_basis(x=domain[0])
                ed_et = self.basises[1].edge_basis(x=domain[1])
                bf_edge_det = np.kron(ed_et, lb_xi)
                ed_xi = self.basises[0].edge_basis(x=domain[0])
                lb_et = self.basises[1].node_basis(x=domain[1])
                bf_edge_dxi = np.kron(lb_et, ed_xi)
                _basis_ = (bf_edge_det, bf_edge_dxi)
            else:
                raise Exception('<2dCSCG Space> : orientation={} is wrong'.format(orientation))
        # ... 2 ...
        elif k == 2:
            bf_xi = self.basises[0].edge_basis(x=domain[0])
            bf_et = self.basises[1].edge_basis(x=domain[1])
            bf = np.kron(bf_et, bf_xi)
            _basis_ = (bf,)
        else:
            raise DimensionError()

        if compute_xietasigma:
            _xietasigma_ = self.___evaluate_basis_functions_meshgrid_2d___(domain)
            return _xietasigma_, _basis_
        else:
            return None, _basis_

    @staticmethod
    def ___evaluate_basis_functions_meshgrid_2d___(domain):
        xi, eta = np.meshgrid(*domain, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        return xi, eta