

from screws.freeze.main import FrozenOnly
from importlib import import_module






class _3dCSCG_VectorField_DO(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._vf_.___DO_evaluate_func_at_time___(time=time)

    def inner_product_with_space_of(self, other, quad_degree=None):
        return self._vf_.___PRIVATE_do_inner_product_with_space_of___(other, quad_degree=quad_degree)

    def reconstruct(self, *args, **kwargs):
        return self._vf_.reconstruct(*args, **kwargs)

    def cross_product(self, vector):
        """We compute self cross_product another vector: self X vector"""

        if vector.__class__.__name__ == '_3dCSCG_VectorField':

            if self._vf_.ftype == 'standard':
                if vector.ftype == 'standard':
                    w0, w1, w2 = self._vf_.func
                    u0, u1, u2 = vector.func
                    CP0 = ___VF_CROSS_PRODUCT_HELPER_1___(w1, u2, w2, u1)
                    CP1 = ___VF_CROSS_PRODUCT_HELPER_1___(w2, u0, w0, u2)
                    CP2 = ___VF_CROSS_PRODUCT_HELPER_1___(w0, u1, w1, u0)

                    vector_class = getattr(import_module('_3dCSCG.fields.vector.main'), '_3dCSCG_VectorField')

                    cp_vector = vector_class(self._vf_.mesh,
                                                    [CP0, CP1, CP2],
                                                    ftype='standard',
                                                    valid_time=self._vf_.valid_time,
                                                    name = self._vf_.standard_properties.name
                                                         + '--cross-X-product--'
                                                         + vector.standard_properties.name
                                                    )
                    return cp_vector
                else:
                    raise NotImplementedError(
                        f"a standard _3dCSCG_VectorField cannot do cross product with a _3dCSCG_VectorField of ftype {vector.ftype}.")
            else:
                raise NotImplementedError(f"a _3dCSCG_VectorField of ftype {self._vf_.ftype} cannot do cross product.")
        else:
            raise NotImplementedError(f"a _3dCSCG_VectorField can not cross product a {vector}.")





class ___VF_CROSS_PRODUCT_HELPER_1___(object):
    def __init__(self, f0, f1, f2, f3):
        self._f0_ = f0
        self._f1_ = f1
        self._f2_ = f2
        self._f3_ = f3

    def __call__(self, t, x, y, z):
        return self._f0_(t, x, y, z) * self._f1_(t, x, y, z) - self._f2_(t, x, y, z) * self._f3_(t, x, y, z)
