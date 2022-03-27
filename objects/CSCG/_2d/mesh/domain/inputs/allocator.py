

import sys
if './' not in sys.path: sys.path.append('/')

from screws.freeze.main import FrozenOnly
from importlib import import_module


class DomainInputFinder(FrozenOnly):
    """ We use this finder to get a `DomainInput`."""

    def __init__(self, ID):
        """
        Parameters
        ----------
        ID : str

        """
        assert ID in self.___defined_DI___(), f" <DomainInputFinder> : mesh ID = {ID} is wrong."
        cls_name = self.___defined_DI___()[ID]
        cls_path = self.___DI_path___()[ID]
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
                  'rectangle': "Rectangle",
                  'rectangle_periodic': "RectanglePeriodic",}
        return _dict_

    @classmethod
    def ___DI_path___(cls):
        """ """
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'

        return {'chp1' : base_path + "chp1",
                'chp2' : base_path + "chp2",
                'crazy': base_path + "crazy",
                'quadrangle'    : base_path + 'quadrangle',
                'crazy_periodic': base_path + "crazy_periodic",
                'bcr': base_path + "bcr",
                'cic': base_path + "cic",
                'rectangle': base_path + "rectangle",
                'rectangle_periodic': base_path + "rectangle_periodic",}







if __name__ == "__main__":
    # mpiexec python _2dCSCG\main.py
    di = DomainInputFinder('rectangle_periodic')()