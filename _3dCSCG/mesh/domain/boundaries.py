"""
This is a class for the boundaries of the domain. A more useful class is the one for the boundaries of the mesh. See
module "boundaries" for the mesh class.

If two boundaries are periodic, in other words, they are internal, we will still have them contained in this domain
boundaries class.

In the mesh.boundaries class, the periodic boundaries will not be shown.
"""

import sys
if './' not in sys.path: sys.path.append('./')

from root.config import *
from SCREWS.frozen import FrozenOnly


class _3dCSCG_Boundaries(FrozenOnly):
    """"""
    def __init__(self, domain):
        assert domain.ndim == 3, " <Domain> <Boundaries> "
        assert domain.__class__.__name__ == '_3dCSCG_Domain', " <Domain> <Boundaries> "
        self._domain_ = domain
        self._distribution_regularities_ = None
        self._freeze_self_()

    @property
    def names(self):
        return self._domain_._boundary_names_

    @property
    def num(self):
        return self._domain_._num_boundaries_

    @property
    def distribution_regularities(self):
        """How the boundaries are distributed. Return a list containing one or some of:

            (1) "Regular:one-region-corner-interface": Any two of boundaries,
                if they connect with each other, should approach each other from
                the same region. Topologically* (regions like orthogonal
                structured mesh cells), all connected boundaries are
                perpendicular (90 degree, NOT 270 degree!!!!) to each other.

        """
        if self._distribution_regularities_ is not None:
            return self._distribution_regularities_

        self._distribution_regularities_ = list()

        if self.___PRIVATE_if_is_Regular__one_region_corner_interface___():
            self._distribution_regularities_.append("Regular:one-region-corner-interface")


        return self._distribution_regularities_


    def ___PRIVATE_if_is_Regular__one_region_corner_interface___(self):
        """"""
        if rAnk == mAster_rank:
            ToF = True
            # for test reasons we
            NUM = self._domain_.regions.num
            if NUM == 1: # only one region, then must be a Regular:one-region-corner-interface
                pass
            else:
                MAP = self._domain_.regions.map

                for rn in MAP:
                    for i, s in enumerate('NSWEBF'): # go through all region sides
                        what_here = MAP[rn][i]
                        if what_here in self.names:
                            if s in 'NS':
                                check_directions = [2,3,4,5]
                            elif s in 'WE':
                                check_directions = [0,1,4,5]
                            elif s in 'BF':
                                check_directions = [0,1,2,3]
                            else:
                                raise Exception()

                            for cd in check_directions:
                                what_is_at_this_side = MAP[rn][cd]
                                if what_is_at_this_side in self.names:
                                    pass # OK, remain True
                                else:
                                    assert what_is_at_this_side[:2] == 'R:', "trivial check."

                                    what_side_connected = MAP[what_is_at_this_side][i]
                                    if what_side_connected == what_here:
                                        pass
                                    else:
                                        ToF = False # only this is happening, we stop.
                                        break

                        else: # internal region side, skip it.
                            assert what_here[:2] == 'R:', "trivial check."

                        if ToF is False: break
                    if ToF is False: break
        else:
            ToF = None
        ToF = cOmm.bcast(ToF, root=mAster_rank)
        return ToF



if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\mesh\domain\boundaries.py
    from _3dCSCG.main import MeshGenerator
    elements = [3, 4, 2]
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    mesh = MeshGenerator('crazy')(elements)

    D = mesh.domain.boundaries.distribution_regularities

    print(D)