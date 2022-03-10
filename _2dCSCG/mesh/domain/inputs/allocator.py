

from screws.freeze.main import FrozenOnly
from screws.exceptions import MeshError
from importlib import import_module


class DomainInputFinder(FrozenOnly):
    """ We use this finder to get a `DomainInput`."""

    def __init__(self, ID):
        """
        Parameters
        ----------
        ID : str

        """
        try:
            mesh_class = self.___defined_DI___()[ID]
        except KeyError:
            raise MeshError(" <DomainInputFinder> : mesh ID = {} is wrong.".format(ID))
        cls_name = mesh_class
        cls_path = self.___DI_path___() + '.' + ID
        self._DomainInput_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """ """
        return self._DomainInput_(*args, **kwargs)

    @classmethod
    def ___defined_DI___(cls):
        """
        Here we store all defined meshComponents. Whenever we define a new meshComponents (
        actually, a new domain_input), we add a nickname for it here.

        """
        _dict_ = {'chp1': "CircleHolePlate1",
                  'chp2': "CircleHolePlate2",
                  'crazy': "Crazy",
                  'quadrangle': 'Quadrangle',
                  'crazy_periodic': "CrazyPeriodic",
                  'bcr': "BottomCustomizedRectangle",
                  'cic': "CylinderInChannel",
                  'rectangle': "Rectangle"}
        return _dict_

    @classmethod
    def ___DI_path___(cls):
        """ """
        return '_2dCSCG.mesh.domain.inputs'


if __name__ == "__main__":
    di = DomainInputFinder('rectangle')()