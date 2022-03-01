
from screws.frozen import FrozenOnly
from root.config import *

class _3dCSCG_Standard_Form_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def reset_cache(self):
        self._sf_.___PRIVATE_reset_cache___()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._sf_.___PRIVATE_do_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def evaluate_basis_at_quadrature(self, quad_degree, quad_type=None, compute_xietasigma=True):
        return self._sf_.space.do.evaluate_form_basis_at_quadrature(
            self._sf_.k, quad_degree, quad_type=quad_type, compute_xietasigma=compute_xietasigma)

    def resemble(self, *args, **kwargs):
        return self._sf_.___PRIVATE_do_resemble___(*args, **kwargs)

    def interpolate(self, data_set):
        """Like resemble, but `interpolate` is from a data set, not another form.

        :param data_set:
        :return:
        """
        raise NotImplementedError()

    def cross_product(self, *args, **kwargs):
        return self._sf_.special.cross_product(*args, **kwargs)

    def compute_L2_inner_product_energy_with(self, *args, **kwargs):
        return self._sf_.___PRIVATE_do_compute_L2_inner_product_energy_with___(*args, **kwargs)

    def compute_L2_diff_from(self, other, n=2, quad_degree=None):
        sf = self._sf_
        of = other
        assert '3dCSCG_standard_form' in of.standard_properties.tags, "Other should be a _3dCSCG standard form."
        if sf.k in (0, 3):
            assert sf.k in (0, 3)
            OneOrThree = 1
        else:
            assert of.k in (1, 2)
            OneOrThree = 3
        if sf.mesh == of.mesh:
            if quad_degree is None: quad_degree = [np.max([sf.dqp[i], of.dqp[i]]) + 1 for i in range(3)]
            quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
            Jacobian = sf.mesh.elements.coordinate_transformation.QUAD_1d.Jacobian(quad_degree, 'Gauss')
            _, SV = sf.reconstruct(*quad_nodes, ravel=True)
            _, OV = of.reconstruct(*quad_nodes, ravel=True)
            localError = list()
            for i in sf.mesh.elements.indices:
                detJ = Jacobian[i]
                LEIntermediate = np.sum([(SV[i][m] - OV[i][m])**n for m in range(OneOrThree)], axis=0)
                localError.append(np.sum(LEIntermediate * detJ * quad_weights))
            core_local = np.sum(localError)
            core_local = cOmm.gather(core_local, root=mAster_rank)
            if rAnk == mAster_rank:
                globalError = np.sum(core_local) ** (1/n)
            else:
                globalError = None
            globalError = cOmm.bcast(globalError, root=mAster_rank)
            return globalError
        else:
            raise NotImplementedError('Can only work on forms from the same mesh now.')

    def discretize(self, *args, **kwargs):
        return self._sf_.discretize(*args, **kwargs)

    def reconstruct(self, *args, **kwargs):
        return self._sf_.reconstruct(*args, **kwargs)
