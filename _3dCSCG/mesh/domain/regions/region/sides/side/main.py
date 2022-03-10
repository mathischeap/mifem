
from screws.freeze.main import FrozenOnly
import numpy as np
from root.config.main import saFe_mode

from _3dCSCG.mesh.domain.regions.region.sides.side.coordinate_transformation import SideCoordinateTransformation

class Side(FrozenOnly):
    def __init__(self, sides, side_name):
        assert sides.ndim == 3, " <Region> <DIII> <Side> "
        assert sides.__class__.__name__ == 'Sides', " <Region> <DIII> <Side> "
        assert side_name in ('N', 'S', 'W', 'E', 'B', 'F'), " <Region> <DIII> <Side> "
        self._sides_ = sides
        self._region_ = sides._region_
        self._position_ = sides._region_.name + '-' + side_name
        self._name_ = side_name
        self._ct_ = None
        self._freeze_self_()

    @property
    def ndim(self):
        return 3

    @property
    def position(self):
        return self._position_



    @property
    def coordinate_transformation(self):
        if self._ct_ is None:
            self._ct_ = SideCoordinateTransformation(self)
        return self._ct_


    def ___generate_full_ep___(self, evaluation_points):
        """ """
        if len(evaluation_points) == 2: # only two valid ep for local trace element provided.

            ep0, ep1 = evaluation_points
            shape0 = np.shape(ep0)
            if saFe_mode:
                assert shape0 == np.shape(ep1), \
                    " <TraceElement3D> : evaluation_points shape wrong."

            if self._name_ == 'N': ep = (np.zeros(shape0), *evaluation_points)
            elif self._name_ == 'S': ep = (np.ones(shape0), *evaluation_points)
            elif self._name_ == 'W': ep = (ep0, np.zeros(shape0), ep1)
            elif self._name_ == 'E': ep = (ep0, np.ones(shape0), ep1)
            elif self._name_ == 'B': ep = (*evaluation_points, np.zeros(shape0))
            elif self._name_ == 'F': ep = (*evaluation_points, np.ones(shape0))
            else:
                raise Exception()

            return ep

        elif len(evaluation_points) == 3: # two valid plus one -1 or +1 coordinates are provided
            if saFe_mode:
                assert evaluation_points[0].shape == \
                       evaluation_points[1].shape == \
                       evaluation_points[2].shape, "evaluation_points shape wrong."
            return evaluation_points

        else:
            raise Exception("evaluation_points shape wrong dimension wrong.")



