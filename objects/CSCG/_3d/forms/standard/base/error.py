# -*- coding: utf-8 -*-

from root.config.main import *
from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.ADF.allocator import _3dCSCG_ADF_Allocator
from importlib import import_module



class _3dCSCG_Standard_Form_Error(FrozenOnly):
    """"""

    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def L(self, n=2, quad_degree=None, upon=False, quad_density=None):
        """The global :math:`L^2` error;

        It is global, so slaves first send info to the secretary who computes the
        global error and sends it back to all slaves.

        :param int n: (`default`:``2``) :math:`L^{n}` error.
        :param quad_degree: The quadrature degree used to compute the error.
        :param bool upon: If True, we will shift all the reconstructed value in order to one of the value
            match the exact value. This is very useful to measure the error like the pressure or potential
            which is not determined before-hands. We have to fix it at one point. We can do it before solving,
            but it is very complicated. So we'd better do ``upon=True`` when measure the error.
        :param quad_density: Only used for n == 'infinity'
        :return: The global :math:`L^{n}` error.
        :rtype: float
        """

        assert self._sf_.CF.ftype == 'standard', \
            f"Currently, this L^n error method only works for standard functions."

        assert self._sf_.cochain.local is not None, " I have no cochain."
        OneOrThree = 1 if self._sf_.k in (0, 3) else 3

        quad_degree = [self._sf_.dqp[i] + 2 for i in range(3)] \
            if quad_degree is None else quad_degree

        if n == 'infinity':
            if quad_density is not None:
                assert isinstance(quad_density, (int, float)) and quad_density > 0, \
                    f"quad_density ={quad_density} must be int or float > 0."
                NUM_elements = self._sf_.mesh.elements.GLOBAL_num
                density_per_element = quad_density / NUM_elements
                num_nodes = density_per_element**(1/3)
                if num_nodes < 1: num_nodes = 3

                if num_nodes % 1 >= 0.5:
                    num_nodes = int(num_nodes) + 1
                else:
                    num_nodes = int(num_nodes)

                _nodes_ = np.linspace(-1, 1, num_nodes+1)
                _nodes_ = (_nodes_[:-1] + _nodes_[1:]) / 2
                quad_nodes = [_nodes_, _nodes_, _nodes_]

            else:
                quad_nodes = [np.linspace(-1, 1, p+10) for p in quad_degree]

        else:
            assert isinstance(n, int) and n > 0, f"L^{n} error is not valid."
            quad_nodes, _, quad_weights = self._sf_.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)



        xi, eta, sigma = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        sigma = sigma.ravel('F')

        xyz, v = self._sf_.reconstruct(*quad_nodes, ravel=True)

        # upon shift ... --------------------------- BELOW --------------------------------------
        add_to = None
        if upon is True: # we will shift according to the first quadrature point.
            if len(xyz) > 0:
                bi = self._sf_.mesh.elements.indices[0]
                bx, by, bz = xyz[bi]
                bv = v[bi]
                base_xyz = [bx[0], by[0], bz[0]]
                if OneOrThree == 1:
                    base_val = [bv[0][0],]
                else:
                    base_val = [bv[0][0], bv[1][0], bv[2][0]]
            else:
                base_xyz = None
                base_val = None
            base_xyz = COMM.gather(base_xyz, root=MASTER_RANK)
            base_val = COMM.gather(base_val, root=MASTER_RANK)
            if RANK == MASTER_RANK:
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
            add_to = COMM.bcast(add_to, root=MASTER_RANK)
        elif upon is False:
            pass
        else:
            raise Exception(f"upon can not be {upon}")
        # ...
        if add_to is not None:
            for i in self._sf_.mesh.elements.indices:
                for j in range(OneOrThree):
                    v[i][j] += add_to[j]
        # +++++++++++++++++++++++++++++++++++++++++++ ABOVE +++++++++++++++++++++++++++++++++++++

        if n == 'infinity':
            localError = -1

            for i in self._sf_.mesh.elements.indices:

                error_i = [np.max(np.abs(v[i][m] -
                    self._sf_.CF.___DO_evaluate_func_at_time___()[m](*xyz[i])))
                           for m in range(OneOrThree)]
                error_i = max(error_i)

                localError = error_i if error_i > localError else localError

            LOC_ERR = COMM.gather(localError, root=MASTER_RANK)

            if RANK == MASTER_RANK:
                globalError = max(LOC_ERR)
            else:
                globalError = None

            globalError = COMM.bcast(globalError, root=MASTER_RANK)

        #--------------- n < infty, n-error------------------------------------------------
        else:

            assert isinstance(n, int) and n > 0, f"L^{n} error is not valid."

            localError = list()
            for i in self._sf_.mesh.elements.indices:
                element = self._sf_.mesh.elements[i]
                detJ = element.coordinate_transformation.Jacobian(xi, eta, sigma)
                LEIntermediate = np.sum(
                    [(v[i][m] - self._sf_.CF.___DO_evaluate_func_at_time___()[m](*xyz[i]))**n
                     for m in range(OneOrThree)], axis=0
                )
                # noinspection PyUnboundLocalVariable
                localError.append(np.sum(LEIntermediate * detJ * quad_weights))

            core_local = np.sum(localError)
            core_local = COMM.gather(core_local, root=MASTER_RANK)

            if RANK == MASTER_RANK:
                globalError = np.sum(core_local) ** (1 / n)
            else:
                globalError = None
            globalError = COMM.bcast(globalError, root=MASTER_RANK)

        #==========================================================================================

        assert globalError >= 0, f"L_{n}error = {globalError} wrong, it must >= 0."

        return globalError

    def H(self, d_func, quad_degree=None, upon=False, quad_density=None):
        """
        Global :math:`H^1`-error; :math:`H^1` error includes :math:`H(\\mathrm{curl})` and :math:`H(\\mathrm{div})`;
        since it basically use the ``globalL2`` method, so the computation is done in the secretary core and
        communications are needed.

        :param d_func: The function of the derivative of ``self`` form.
        :param quad_degree: The quadrature degree used to compute the error.
        :param upon:
        :param quad_density:
        :return: The global :math:`H^{1}` error.
        :rtype: float
        """
        selfErrorL2 = self.L(n=2, quad_degree=quad_degree, upon=upon, quad_density=quad_density)
        D_self = self._sf_.coboundary()
        D_self.CF = d_func
        D_self.CF.current_time = self._sf_.CF.current_time
        DErrorL2 = D_self.error.L(n=2, quad_degree=quad_degree)
        return (selfErrorL2 ** 2 + DErrorL2 ** 2) ** 0.5

    def dH(self, dt, dfunc, time=None, n=1):
        """"""

        k = self._sf_.k

        dual_class_name = _3dCSCG_ADF_Allocator.___forms_name___()[str(k) + '-adf']
        dual_class_path = _3dCSCG_ADF_Allocator.___forms_path___()[str(k) + '-adf']

        CLASS = getattr(import_module(dual_class_path), dual_class_name)

        sf = self._sf_
        mesh = sf.mesh
        space = sf.space
        orientation = sf.orientation
        name = 'dual-' +  sf.standard_properties.name

        if time is None:
            time = self._sf_.CF.current_time

        dual_form = CLASS(self._sf_, mesh, space, orientation=orientation, name=name)

        return dual_form.error.dH(dt, dfunc, time=time, n=n)