"""
Sub-geometries of a regions.

For example, we can pick a slice from a regions. Then we can do reconstruction of a form on this slice to view some
particular structure.

"""

import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.domain.regions.region.sub_geometry.perpendicular_slice import RegionPerpendicularSlice


class RegionSubGeometry(FrozenOnly):
    """"""

    def __init__(self, region):
        """"""
        self._region_ = region

        self._freeze_self_()





    def make_a_perpendicular_slice_object_on(self, x=None, y=None, z=None, r=None, s=None, t=None):
        """"""
        num_None = 0
        for _ in (x, y, z, r, s, t):
            if _ is None:
                num_None += 1
        assert num_None == 5, f"please only provide one of x, y, z, r, s, t."

        if any([_ is not None for _ in (x, y, z)]):
            if x is not None:
                assert isinstance(x, (int, float)), f"I need a float or int"
                r = self._region_.interpolation.___inverse_mapping_r_x_s0t0___(x)
                r = float(r)
            if y is not None:
                assert isinstance(y, (int, float)), f"I need a float or int"
                s = self._region_.interpolation.___inverse_mapping_s_y_r0t0___(y)
                s = float(s)
            if z is not None:
                assert isinstance(z, (int, float)), f"I need a float or int"
                t = self._region_.interpolation.___inverse_mapping_t_z_r0s0___(z)
                t = float(t)

            for _ in (r, s, t):
                if _ is not None:
                    if 0 <= _ <= 1:
                        pass
                    else:
                        return None # this region has no business with this slice x, y or z constant slice.

            return RegionPerpendicularSlice(self._region_, r=r, s=s, t=t)

        elif any([_ is not None for _ in (r, s, t)]):
            if r is not None: assert 0<= r <= 1
            if s is not None: assert 0<= s <= 1
            if t is not None: assert 0<= t <= 1
            return RegionPerpendicularSlice(self._region_, r=r, s=s, t=t)

        else:
            raise Exception()




if __name__ == '__main__':
    # mpiexec -n 8 python objects/CSCG/_3d/mesh/domain/regions/region/sub_geometry/main.py

    from objects.CSCG._3d.master import MeshGenerator

    mesh = MeshGenerator('crazy', bounds=([0,2],[0,2],[0,2]))([3, 3, 3], EDM='chaotic', show_info=True)

    R = mesh.domain.regions['R:R']


    RSG = R.sub_geometry
    RSG.make_a_perpendicular_slice_object_on(x=2)

