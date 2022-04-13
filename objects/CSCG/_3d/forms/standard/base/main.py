# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from objects.CSCG._3d.forms.base import _3dCSCG_FORM_BASE
from objects.CSCG._3d.forms.standard.base.numbering.main import _3dCSCG_Standard_Form_Numbering
from objects.CSCG._3d.forms.standard.base.export.main import _3dCSC_Standard_Form_Export
from root.config.main import *
from root.save import read
from scipy.interpolate import NearestNDInterpolator

from objects.CSCG.base.forms.standard.main import CSCG_Standard_Form

from objects.CSCG._3d.forms.standard.base.dofs.main import _3dCSCG_Standard_forms_DOFs
from objects.CSCG._3d.forms.standard.base.operators.main import _3dCSCG_Standard_Form_Operators
from objects.CSCG._3d.forms.standard.base.matrices import _3dCSCG_Standard_Form_Matrices
from objects.CSCG._3d.forms.standard.base.coboundary import _3dCSCG_Standard_Form_Coboundary
from objects.CSCG._3d.forms.standard.base.error import _3dCSCG_Standard_Form_Error
from objects.CSCG._3d.forms.standard.base.cochain.main import _3dCSCG_Standard_Form_Cochain
from objects.CSCG._3d.forms.standard.base.do import _3dCSCG_Standard_Form_DO



class _3dCSCG_Standard_Form(CSCG_Standard_Form, _3dCSCG_FORM_BASE, ndim=3):
    """
    This is the parent of all 3d standard forms.

    :param mesh:
    :param space:
    :param bool is_hybrid:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as the scheme name, and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, is_hybrid, orientation, numbering_parameters, name):
        super().__init__(mesh, space)
        super().___init___()
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
        self._DO_ = _3dCSCG_Standard_Form_DO(self)

        self._export_ = _3dCSC_Standard_Form_Export(self)
        self._dofs_ = None


    def ___PRIVATE_reset_cache___(self):
        self.cochain.___PRIVATE_reset_cache___()
        self.coboundary.___PRIVATE_reset_cache___()

    @property
    def numbering(self):
        return self._numbering_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def error(self):
        return self._error_

    @property
    def coboundary(self):
        return self._coboundary_

    @property
    def matrices(self):
        return self._matrices_

    @property
    def operators(self):
        return self._operators_

    @property
    def do(self):
        """If it has too many do methods, we group them in to do."""
        return self._DO_

    @property
    def export(self):
        return self._export_

    @property
    def dofs(self):
        if self._dofs_ is None:
            self._dofs_ = _3dCSCG_Standard_forms_DOFs(self)
        return self._dofs_





    def ___PRIVATE_do_evaluate_basis_at_meshgrid___(self, xi, eta, sigma, compute_xietasigma=True):
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
        return self.space.do.evaluate_form_basis_at_meshgrid(
            self.k, xi, eta, sigma, compute_xietasigma=compute_xietasigma)

    def ___PRIVATE_do_resemble___(self, obj_or_filename, density=80000):
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
            self.discretize._standard_()
            # self.___PRIVATE_discretize_standard_ftype___()
            self.func._body_ = None

    def ___PRIVATE_saving_info___(self):
        """"""
        my_info = dict()
        my_info['name'] = self.__class__.__name__
        my_info['mesh_ID'] = str(self.mesh)
        return my_info