# -*- coding: utf-8 -*-


from objects.CSCG._3d.mesh.domain.base import _3dCSCG_DomainBase
from objects.CSCG._3d.mesh.domain.visualize import _3dCSCG_Domain_Visualize
from objects.CSCG._3d.mesh.domain.boundaries import _3dCSCG_Boundaries
from objects.CSCG._3d.mesh.domain.IS import _3dCSCG_Domain_IS




class _3dCSCG_Domain(_3dCSCG_DomainBase):
    """We have the whole ``_3dCSCG_Domain`` (all same) in all cores. This
    is very important.
    """
    def __init__(self, di):
        super(_3dCSCG_Domain, self).__init__(di)
        v = 0
        for rn in self.regions.names:
            v += self.regions(rn).volume
        self._volume_ = v
        self._visualize_ = None
        self._boundaries_ = None
        self.___define_parameters___ = None
        self._IS_ = _3dCSCG_Domain_IS(self)
        for rn in self.regions:
            self.regions[rn]._MAP_ = self._region_map_[rn]
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
    def domain_input(self):
        return self._domain_input_

    @property
    def ndim(self):
        return self.domain_input.ndim

    @property
    def volume(self):
        return self._volume_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_Domain_Visualize(self)
        return self._visualize_

    @property
    def boundaries(self):
        if self._boundaries_ is None:
            self._boundaries_ = _3dCSCG_Boundaries(self)
        return self._boundaries_

    @property
    def interpolators(self):
        """
        Returns
        -------
        self._interpolator_ : dict

        """
        return self._interpolators_

    @property
    def name(self):
        return self.domain_input.domain_name

    def __len__(self):
        """
        The length of a domain, len(domain), is defined to be the number of
        regions the domain has.

        """
        return self._num_regions_

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