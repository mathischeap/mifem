# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from importlib import import_module
from screws.frozen import FrozenOnly
from _3dCSCG.forms.base import _3dCSCG_FORM_BASE
from _3dCSCG.forms.standard.base.numbering.main import _3dCSCG_Standard_Form_Numbering
from _3dCSCG.forms.standard.base.visualize.main import _3dCSCG_FormVisualize
from _3dCSCG.forms.standard.base.export.main import _3dCSC_Standard_Form_Export
from root.config import *
from tools.linear_algebra.elementwise_cache import EWC_SparseMatrix
from root.mifem import read
from scipy.interpolate import NearestNDInterpolator

from inheriting.CSCG.form.standard.main_BASE import CSCG_Standard_Form
from inheriting.CSCG.form.standard.main_BASE import CSCG_Standard_Form_Cochain_BASE
from inheriting.CSCG.form.standard.main_BASE import CSCG_Standard_Form_Coboundary_BASE

from _3dCSCG.forms.standard.base.dofs.main import _3dCSCG_Standard_forms_DOFs

class _3dCSCG_Standard_Form(CSCG_Standard_Form, _3dCSCG_FORM_BASE, ndim=3):
    """
    This is the parent of all 3d standard forms.

    :param mesh:
    :param space:
    :param bool is_hybrid:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as scheme name and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, is_hybrid, orientation, numbering_parameters, name):
        super().__init__(mesh, space)
        self._NUM_basis_, self._NUM_basis_components_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        assert isinstance(is_hybrid, bool), " isHybrid needs to be bool."
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._IS_hybrid_ = is_hybrid
        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_form')
        self._numbering_ = _3dCSCG_Standard_Form_Numbering(self, numbering_parameters)
        self._cochain_ = _3dCSCG_Standard_Form_Cochain(self)
        self._error_ = _3dCSCG_Standard_Form_Error(self)
        self._coboundary_ = _3dCSCG_Standard_Form_Coboundary(self)


        self._matrices_ = _3dCSCG_Standard_Form_Matrices(self)
        self._operators_ = _3dCSCG_Standard_Form_Operators(self)
        self._visualize_ = _3dCSCG_FormVisualize(self)
        self._export_ = _3dCSC_Standard_Form_Export(self)
        self._DO_ = _3dCSCG_Standard_Form_DO(self)
        self._dofs_ = None

    def RESET_cache(self):
        self.cochain.RESET_cache()
        self.coboundary.RESET_cache()

    def ___DO_evaluate_basis_at_meshgrid___(self, xi, eta, sigma, compute_xietasigma=True):
        """
        Evaluate the basis functions on ``meshgrid(xi, eta, sigma)``.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param sigma: A 1d iterable object of floats between -1 and 1.
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :type sigma: list, tuple, numpy.ndarray
        :param bool compute_xietasigma: (`default`:``True``) If we compute the
            ``meshgrid(xi, eta, sigma, indexing='ij')``.
        :returns: A tuple of outputs:

            1. (None, tuple) -- ``(xi, eta, sigma)`` after ``meshgrid`` and ``ravel('F')``.
            2. tuple -- The evaluated basis functions.
        """
        return self.space.DO_evaluate_form_basis_at_meshgrid(
            self.k, xi, eta, sigma, compute_xietasigma=compute_xietasigma)

    def ___DO_resemble___(self, obj_or_filename, density=80000):
        """
        :param obj_or_filename:
        :param density:
        :return:
        """
        if isinstance(obj_or_filename, str):
            of = read(obj_or_filename)
        else:
            of = obj_or_filename

        assert self.mesh.domain == of.mesh.domain, "domain must be same."
        assert self.__class__.__name__ == of.__class__.__name__
        assert self.mesh.__class__.__name__ == of.mesh.__class__.__name__

        if self.mesh == of.mesh and self.space == of.space:
            # this is the simplest case, just copy the cochain.
            self.cochain.local = of.cochain.local
        else:
            bp = int(np.ceil((density / self.mesh.elements.GLOBAL_num) ** (1/3)))
            p = [bp + self.p[i] for i in range(3)]
            gap = [1 / (p[i]+1) for i in range(3)]
            r = np.linspace(-1 + gap[0], 1 - gap[0], p[0])
            s = np.linspace(-1 + gap[1], 1 - gap[1], p[1])
            t = np.linspace(-1 + gap[2], 1 - gap[2], p[2])
            xyz, V = of.reconstruct(r, s, t, ravel=True)
            LEN = 1 if self.k in (0, 3) else 3
            xyz = cOmm.gather(xyz, root=mAster_rank)
            V = cOmm.gather(V, root=mAster_rank)
            if rAnk == mAster_rank:
                XYZ = dict()
                VVV = dict()
                for i in range(len(xyz)):
                    XYZ.update(xyz[i])
                    VVV.update(V[i])
                del xyz, V
                X = list()
                Y = list()
                Z = list()
                V = [list() for _ in range(LEN)]
                for i in range(of.mesh.elements.GLOBAL_num):
                    X.extend(XYZ[i][0])
                    Y.extend(XYZ[i][1])
                    Z.extend(XYZ[i][2])
                    for j in range(LEN):
                        V[j].extend(VVV[i][j])
                del XYZ, VVV
                for i in range(LEN):
                    # noinspection PyTypeChecker
                    V[i] = np.array(V[i])
                X = np.array(X)
                Y = np.array(Y)
                Z = np.array(Z)
                func = list()
                for i in range(LEN):
                    func.append(NearestNDInterpolator((X, Y, Z), V[i]))
            else:
                func = None
            func = cOmm.bcast(func, root=mAster_rank)
            self.func._body_ = func
            # noinspection PyUnresolvedReferences
            self.___PRIVATE_discretize_standard_ftype___()
            self.func._body_ = None

    def ___DO_compute_L2_inner_product_energy_with___(self, other=None, M=None):
        """
        Compute (self, other)_{L2}.

        :param other:
        :param M:
        :return:
        """
        if other is None: other = self
        assert self.mesh == other.mesh, "Meshes do not match."
        if M is None: M = self.operators.inner(other)
        LOCAL = list()
        for i in self.mesh.elements:
            LOCAL.append(self.cochain.local[i] @ M[i] @ other.cochain.local[i])
        LOCAL = np.sum(LOCAL)
        return cOmm.allreduce(LOCAL, op=MPI.SUM)

    @property
    def export(self):
        return self._export_

    @property
    def dofs(self):
        if self._dofs_ is None:
            self._dofs_ = _3dCSCG_Standard_forms_DOFs(self)
        return self._dofs_


    def ___PRIVATE_saving_info___(self):
        """"""
        my_info = dict()
        my_info['name'] = self.__class__.__name__
        my_info['mesh_ID'] = str(self.mesh)
        return my_info


class _3dCSCG_Standard_Form_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def RESET_cache(self):
        self._sf_.RESET_cache()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._sf_.___DO_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def evaluate_basis_at_quadrature(self, quad_degree, quad_type=None, compute_xietasigma=True):
        return self._sf_.space.DO_evaluate_form_basis_at_quadrature(
            self._sf_.k, quad_degree, quad_type=quad_type, compute_xietasigma=compute_xietasigma)

    def resemble(self, *args, **kwargs):
        return self._sf_.___DO_resemble___(*args, **kwargs)

    def interpolate(self, data_set):
        """Like resemble, but interpolate is from a data set, not another form.

        :param data_set:
        :return:
        """
        raise NotImplementedError()

    def cross_product(self, *args, **kwargs):
        return self._sf_.special.cross_product(*args, **kwargs)

    def compute_L2_inner_product_energy_with(self, *args, **kwargs):
        return self._sf_.___DO_compute_L2_inner_product_energy_with___(*args, **kwargs)

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
            quad_nodes, _, quad_weights = sf.space.DO_evaluate_quadrature(quad_degree)
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



class _3dCSCG_Standard_Form_Cochain(CSCG_Standard_Form_Cochain_BASE):
    """The cochain must be full. So it must represent all dofs. For partial dofs, we can access them
    through `cochain.partial` property.
    """
    def __init__(self, sf):
        self._partial_ = None
        super().__init__(sf)

    def local_(self, axis):
        """
        The local cochain along a particular axis.

        :param str axis: The local cochain along which axis? ``x``, ``y`` or ``z``.
        :return: The local cochain dict.
        :rtype: Dict[int, numpy.ndarray]
        """
        assert self._sf_.k not in (0, 3), \
            " <Cochain> : %r is a scalar form; has no idea of axes, use cochain.globe." % self._sf_
        numOfBasisComponents = self._sf_.NUM_basis_components
        localAlongAxis = dict()
        for i in self._sf_.mesh.elements:
            if axis == 'x':
                localAlongAxis[i] = self.local[i][:numOfBasisComponents[0]]
            elif axis == 'y':
                localAlongAxis[i] = self.local[i][numOfBasisComponents[0]:
                                                  numOfBasisComponents[0]+numOfBasisComponents[1]]
            elif axis == 'z':
                localAlongAxis[i] = self.local[i][-numOfBasisComponents[2]:]
            else:
                raise Exception()
        return localAlongAxis


    @property
    def partial(self):
        """To access partial of the cochain."""
        if self._partial_ is None:
            self._partial_ = _3dCSCG_Standard_Form_Cochain_Partial(self)
        return self._partial_


class _3dCSCG_Standard_Form_Cochain_Partial(FrozenOnly):
    """For accessing a partial of the cochain."""
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()










class _3dCSCG_Standard_Form_Error(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def L(self, n=2, quad_degree=None, upon=False, quad_density=None):

        """
        The global :math:`L^2` error; it is global, so slaves first send info to the secretary who computes the
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

        assert self._sf_.func.ftype == 'standard', f"Currently, this L^n error method only works for standard functions."

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
            quad_nodes, _, quad_weights = self._sf_.space.DO_evaluate_quadrature(quad_degree)



        xi, eta, sigma = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        sigma = sigma.ravel('F')

        xyz, v = self._sf_.reconstruct(*quad_nodes, ravel=True)

        # upon shift ... --------------------------- BELOW -------------------------------------------------------------
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
        # +++++++++++++++++++++++++++++++++++++++++++ ABOVE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if n == 'infinity':
            localError = -1

            for i in self._sf_.mesh.elements.indices:

                error_i = [np.max(np.abs(v[i][m] - self._sf_.func.body[m](*xyz[i]))) for m in range(OneOrThree)]
                error_i = max(error_i)

                localError = error_i if error_i > localError else localError

            LOC_ERR = cOmm.gather(localError, root=mAster_rank)

            if rAnk == mAster_rank:
                globalError = max(LOC_ERR)
            else:
                globalError = None

            globalError = cOmm.bcast(globalError, root=mAster_rank)

        else:

            assert isinstance(n, int) and n > 0, f"L^{n} error is not valid."

            localError = list()
            for i in self._sf_.mesh.elements.indices:
                element = self._sf_.mesh.elements[i]
                detJ = element.coordinate_transformation.Jacobian(xi, eta, sigma)
                LEIntermediate = np.sum(
                    [(v[i][m] - self._sf_.func.body[m](*xyz[i]))**n for m in range(OneOrThree)], axis=0
                )
                # noinspection PyUnboundLocalVariable
                localError.append(np.sum(LEIntermediate * detJ * quad_weights))

            core_local = np.sum(localError)
            core_local = cOmm.gather(core_local, root=mAster_rank)

            if rAnk == mAster_rank:
                globalError = np.sum(core_local) ** (1 / n)
            else:
                globalError = None
            globalError = cOmm.bcast(globalError, root=mAster_rank)

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
        D_self.TW.func.body = d_func
        D_self.TW.current_time = self._sf_.TW.current_time
        D_self.TW.DO.push_all_to_instant()
        DErrorL2 = D_self.error.L(n=2, quad_degree=quad_degree)
        return (selfErrorL2 ** 2 + DErrorL2 ** 2) ** 0.5



class _3dCSCG_Standard_Form_Coboundary(CSCG_Standard_Form_Coboundary_BASE):
    def __init__(self, sf):
        super().__init__(sf)

    @property
    def incidence_matrix(self):
        """(scipy.sparse.csc_matrix) Return ths incidence matrix of the standard form."""
        assert self._sf_.k < 3, "volume form has no incidence matrix."
        if self._incidenceMatrix_ is None:
            formName = self._sf_.__class__.__name__
            _incidenceMatrix_ = getattr(self._sf_.space.incidence_matrix, formName)
            self._incidenceMatrix_ = \
                EWC_SparseMatrix(self._sf_.mesh.elements, _incidenceMatrix_, 'constant')

        return self._incidenceMatrix_

    def ___PRIVATE_next_class___(self):
        assert self._sf_.k < 3, "volume form has no next prime space."
        k = self._sf_.k
        nextPath = f'_3dCSCG.forms.standard._{k+1}_form'
        nextName = f'_3dCSCG_{k+1}Form'
        return getattr(import_module(nextPath), nextName)



class _3dCSCG_Standard_Form_Matrices(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._T_ = None
        self._S_ = None
        self._freeze_self_()

    @property
    def mass(self):
        """(Dict[int, scipy.sparse.csr_matrix]) The mass matrix."""
        M = self._sf_.operators.inner(self._sf_)
        # note that even all mesh elements are unique, we still cache the local mass matrices because we may use them for multiple times.
        M.gathering_matrices = (self._sf_, self._sf_)
        return M

    @property
    def incidence(self):
        return self._sf_.coboundary.incidence_matrix

    @property
    def trace(self):
        """Return the trace matrix."""
        if self._T_ is None:
            k = self._sf_.k
            formName = f'_3dCSCG_{int(k)}Trace'
            T = getattr(self._sf_.space.trace_matrix, formName)[0]
            self._T_ = \
                EWC_SparseMatrix(self._sf_.mesh.elements, T, 'constant')
        return self._T_


class _3dCSCG_Standard_Form_Operators(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def inner(self, other, quad_degree=None):
        """
        We do ``(self, other)``, and note that here we only return a matrix; we do not do the inner product
        which needs that both forms have cochain.

        :param other: The other form.
        :param quad_degree:
        """
        data_generator = ___Operators_Inner___(self._sf_, other, quad_degree=quad_degree)
        # note that even all mesh elements are unique, we still cache the output because we may use it for multiple times.
        return EWC_SparseMatrix(self._sf_.mesh.elements, data_generator,
                                self._sf_.mesh.elements.___PRIVATE_elementwise_cache_metric_key___)

    def wedge(self, other, quad_degree=None):
        """"""
        data = self._sf_.___OPERATORS_wedge___(other, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, data, 'constant')

    def cross_product(self, *args, **kwargs):
        return self._sf_.special.cross_product(*args, **kwargs)





class ___Operators_Inner___(FrozenOnly):
    """The class for the inner product matrix."""
    def __init__(self, sf, of, quad_degree=None):
        assert sf.ndim == of.ndim and sf.k == of.k, " <___STORAGE_OPERATORS_INNER___> "
        assert sf.mesh == of.mesh, "Meshes do not match."
        self._mesh_ = sf.mesh
        self._sf_ = sf
        self._of_ = of
        if quad_degree is None:
            quad_degree = [int(np.max([sf.dqp[i], of.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = sf.space.DO_evaluate_quadrature(quad_degree)
        xietasigma, bfSelf = sf.DO.evaluate_basis_at_meshgrid(*quad_nodes)
        if of is sf:
            bfOther = bfSelf
        else:
            xietasigma, bfOther = of.DO.evaluate_basis_at_meshgrid(*quad_nodes)
        self._xietasigma_ = xietasigma
        self._quad_weights_ = quad_weights
        self._bfSelf_ = bfSelf
        self._bfOther_ = bfOther
        self._freeze_self_()

    def __call__(self, i):
        Mi = self._sf_.___OPERATORS_inner___(
            self._of_, i, self._xietasigma_, self._quad_weights_, self._bfSelf_, self._bfOther_
        )
        return Mi
