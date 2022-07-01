# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly

from root.config.main import np, cOmm, rAnk, mAster_rank


class _2dCSCG_Standard_Form_Error(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def L(self, n=2, quad_degree=None, upon=False):
        """
        The global :math:`L^2` error; it is global, so slaves first send info to the secretary who computes the
        global error and sends it back to all slaves.

        :param int n: (`default`:``2``) :math:`L^{n}` error.
        :param quad_degree: The quadrature degree used to compute the error.
        :param bool upon: If True, we will shift all the reconstructed value in order to one of the value
            match the exact value. This is very useful to measure the error like the pressure or potential
            which is not determined before-hands. We have to fix it at one point. We can do it before solving,
            but it is very complicated. So we'd better do ``upon=True`` when measure the error.
        :return: The global :math:`L^{n}` error.
        :rtype: float
        """
        assert self._sf_.cochain.local is not None, " I have no cochain."
        OneOrThree = 1 if self._sf_.k in (0, 2) else 2
        quad_degree = [self._sf_.dqp[i] + 2 for i in range(2)] \
            if quad_degree is None else quad_degree
        quad_nodes, _, quad_weights = self._sf_.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        xi, eta = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        xyz, v = self._sf_.reconstruct(*quad_nodes, ravel=True)

        # upon shift ...
        add_to = None
        if upon is True: # we will shift according to the first quadrature point.
            if len(xyz) > 0:
                bi = self._sf_.mesh.elements.indices[0]
                bx, by = xyz[bi]
                bv = v[bi]
                base_xyz = [bx[0], by[0]]
                if OneOrThree == 1:
                    base_val = [bv[0][0],]
                else:
                    base_val = [bv[0][0], bv[1][0]]
            else:
                base_xyz = None
                base_val = None
            base_xyz = cOmm.gather(base_xyz, root=mAster_rank)
            base_val = cOmm.gather(base_val, root=mAster_rank)
            if rAnk == mAster_rank:
                for i, bi in enumerate(base_xyz):
                    if bi is not None:
                        break
                # noinspection PyUnboundLocalVariable
                assert bi is not None
                # noinspection PyUnboundLocalVariable
                base_xyz = bi # find the base_xyz
                # noinspection PyUnboundLocalVariable
                base_val = base_val[i]
                func_val = [self._sf_.func.body[m](*base_xyz) for m in range(OneOrThree)]
                assert len(base_val) == len(func_val)
                add_to = [func_val[m] - base_val[m] for m in range(OneOrThree)]
                # this add_to will be added to all reconstructed value.
            else:
                pass
            add_to = cOmm.bcast(add_to, root=mAster_rank)
        elif upon is False:
            pass
        else:
            raise Exception(f"upon can not be {upon}")
        # ...
        if add_to is not None:
            for i in self._sf_.mesh.elements.indices:
                for j in range(OneOrThree):
                    v[i][j] += add_to[j]
        # ...

        localError = list()

        for i in self._sf_.mesh.elements.indices:
            element = self._sf_.mesh.elements[i]
            detJ = element.coordinate_transformation.Jacobian(xi, eta)
            LEIntermediate = np.sum(
            [(v[i][m] - self._sf_.func.body[m](*xyz[i]))**n for m in range(OneOrThree)], axis=0
            )
            localError.append(np.sum(LEIntermediate * detJ * quad_weights))

        core_local = np.sum(localError)
        core_local = cOmm.gather(core_local, root=mAster_rank)

        if rAnk == mAster_rank:
            globalError = np.sum(core_local) ** (1 / n)
        else:
            globalError = None
        globalError = cOmm.bcast(globalError, root=mAster_rank)

        return globalError

    def H(self, d_func, quad_degree=None, upon=False):
        """
        Global :math:`H^1`-error; :math:`H^1` error includes :math:`H(\mathrm{curl})` and :math:`H(\mathrm{div})`;
        since it basically use the ``globalL2`` method, so the computation is done in the secretary core and
        communications are needed.

        :param d_func: The function of the derivative of ``self`` form.
        :param quad_degree: The quadrature degree used to compute the error.
        :param upon:
        :return: The global :math:`H^{1}` error.
        :rtype: float
        """
        selfErrorL2 = self.L(n=2, quad_degree=quad_degree, upon=upon)
        D_self = self._sf_.coboundary()
        D_self.TW.func.body = d_func
        D_self.TW.current_time = self._sf_.TW.current_time
        D_self.TW.___DO_push_all_to_instant___()
        DErrorL2 = D_self.error.L(n=2, quad_degree=quad_degree)
        return (selfErrorL2 ** 2 + DErrorL2 ** 2) ** 0.5


