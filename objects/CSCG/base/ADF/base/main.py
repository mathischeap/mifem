# -*- coding: utf-8 -*-
"""
Every dual form is built on a prime form. So, the most important property of a dual form is
"prime" which returns the
prime of it.

"""


from screws.freeze.main import FrozenClass
from root.config.main import *
import pickle

from objects.CSCG.base.ADF.base.BC.main import CSCG_ADForm_BC


class CSCG_Algebra_DUAL_FORM_BASE(FrozenClass):
    """"""
    def __init__(self, ndim, mesh, space, prime):
        """"""
        self.___define_parameters___ = None
        self._prime_ = prime
        self.prime.___is_wrapped_in_ADF___ = True
        self._orientation_ = None
        self.standard_properties.___PRIVATE_add_tag___('CSCG_algebra_dual_form')

        self._mesh_ = mesh
        self._space_ = space

        self.___ndim___ = ndim
        self._k_ = None

        self._mass_matrix_ = None
        self._inverse_mass_matrix_ = None

        self._BC_ = None


    @property
    def ___IS_ADF___(self):
        return True

    @property
    def ___parameters___(self):
        """___parameters___ will be called when you call standard_properties.parameters.

        This is mandatory for FrozenClass.

        :return:
        """
        parameters = dict()
        parameters.update(self.___define_parameters___)
        return parameters

    def ___PRIVATE_save___(self, filename, do_save=False):
        """We have a special saving scheme for AD forms. We basically save its prime form...

        So we write this private save method to override the default save method.

        Better be called from `mifem.save` when save an object (better to save an object by
        using `mifem.save(...)`).
        """
        _2bs_ = dict()
        _2bs_['obj'] = self.__class__.__name__
        _2bs_['parameters'] = self.prime.standard_properties.parameters # this is the key
        if rAnk == mAster_rank:
            if do_save:
                if filename[-3:] != '.mi': filename += '.mi'
                assert filename.count('.') == 1, f"filename={filename} wrong."
                with open(filename, 'wb') as output:
                    pickle.dump(_2bs_, output, pickle.HIGHEST_PROTOCOL)
                output.close()
        return _2bs_

    #-----------------------------------------------------------------------------
    @property
    def BC(self):
        if self._BC_ is None:
            self._BC_ = CSCG_ADForm_BC(self)
        return self._BC_

    @property
    def prime(self):
        """Return the prime form."""
        return self._prime_

    #-----------------------------------------------------------------------------

    @property
    def orientation(self):
        return self._orientation_



    @property
    def mesh(self):
        """Return the mesh."""
        return self._mesh_

    @property
    def space(self):
        """Return the basis function space."""
        return self._space_

    @property
    def p(self):
        """Return the degree of basis functions."""
        return self.space.p



    @property
    def ndim(self):
        """Return the dimensions."""
        return self.___ndim___

    @property
    def k(self):
        return self._k_



    @property
    def mass_matrix(self):
        """For algebra dual forms, we cache the mass matrix. While it is not cached for prime forms."""
        if self._mass_matrix_ is None:
            self._mass_matrix_, self._inverse_mass_matrix_ = self.___PRIVATE_generate_mass_matrix___()
        return self._mass_matrix_

    @property
    def inverse_mass_matrix(self):
        """For algebra dual forms, we cache the mass matrix. While it is not cached for prime forms."""
        if self._inverse_mass_matrix_ is None:
            self._mass_matrix_, self._inverse_mass_matrix_ = self.___PRIVATE_generate_mass_matrix___()
        return self._inverse_mass_matrix_

    def ___PRIVATE_generate_mass_matrix___(self):
        raise NotImplementedError()