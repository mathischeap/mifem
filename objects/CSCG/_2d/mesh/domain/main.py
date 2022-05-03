# -*- coding: utf-8 -*-
"""INTRO

@author: Yi Zhang. Created on Tue May 21 14:10:36 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft,
         Delft, the Netherlands

"""
from objects.CSCG._2d.mesh.domain.base import _2dCSCG_DomainBase
from objects.CSCG._2d.mesh.domain.visualize import _2dCSCG_Domain_Visualize
from objects.CSCG._2d.mesh.domain.boundaries.main import _2dCSCG_Domain_Boundaries
from objects.CSCG._2d.mesh.domain.IS import _2dCSCG_Domain_IS


class _2dCSCG_Domain(_2dCSCG_DomainBase):
    def __init__(self, di):
        """
        Parameters
        ---------
        di : DomainInput
            The DomainInput instance.

        """
        super(_2dCSCG_Domain, self).__init__(di)
        self._visualize_ = _2dCSCG_Domain_Visualize(self) # will only do thing in master core.
        self._boundaries_ = _2dCSCG_Domain_Boundaries(self)
        self.___define_parameters___ = None
        self._IS_ = _2dCSCG_Domain_IS(self)
        self._freeze_self_()

    @property
    def parameters(self):
        return self.___define_parameters___

    def __eq__(self, other):
        for key in self.parameters:
            if key in self.domain_input.internal_parameters:
                pass
            else:
                if self.parameters[key] != other.parameters[key]:
                    return False
        return True


    @property
    def ndim(self):
        return self.domain_input.ndim

    @property
    def domain_input(self):
        return self._domain_input_

    @property
    def name(self):
        return self.domain_input.domain_name

    @property
    def visualize(self):
        return self._visualize_

    @property
    def boundaries(self):
        return self._boundaries_


    def __len__(self):
        """
        The length of a domain, len(domain), is defined to be the number of
        regions the domain has.

        """
        return self._num_regions_


    @property
    def interpolators(self):
        """
        Returns
        -------
        self._interpolator_ : dict

        """
        return self._interpolators_

    @property
    def regions(self):
        """
        Returns
        -------
        self._regions_ : dict
            A dict that contains all regions. Keys: regions' names. Values: the
            regions instances.

        """
        return self._regions_

    @property
    def IS(self):
        return self._IS_