# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/20 11:40 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import np
from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.refinements.future.h_adapt.s0f import hAdapt_S0F


class mpRfT2_Mesh_FutureRefinements_hAdapt(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, obj, approach=1, **kwargs):
        """Let the future mesh adapt to an obj."""
        assert obj.mesh is self._mesh_, f"meshes are different."

        if obj.__class__.__name__ in ('mpRfT2_Si0F', 'mpRfT2_So0F'):
            hAdapt_S0F(obj)(approach=approach, **kwargs)
        else:
            raise NotImplementedError()





if __name__ == "__main__":

    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/refinements/future/h_adapt/main.py

    from __init__ import rfT2
    from objects.mpRfT._2d.forms.standard._0.inner.main import mpRfT2_Si0F
    from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar

    from screws.generators.counter import Counter
    IC = Counter()
    image_folder = '__h_refinement_demo_s0f'
    from screws.miscellaneous.mios import mkdir
    mkdir(image_folder)
    from screws.video.make_a_video_from_images_in_folder import make_a_video_from_images_in_folder

    mesh = rfT2.mesh('crazy', c=0.1, bounds=([-2, 2], [-2, 2]))([15, 15], 3)

    T = np.concatenate([np.linspace(2.5, 0.2, 100), np.linspace(0.2, 2.5, 100)])

    for i in T:
        def p(t, x, y): return 1 / np.exp(np.abs(x**2 + y**2 - i))**2 + 0 * t
        s = mpRfT2_Scalar(mesh, p)
        f = mpRfT2_Si0F(mesh)
        f.TW.func = s
        s.current_time = 0
        f.discretization()
        f.visualization(saveto=image_folder + '/' + str(next(IC)), dpi=150, usetex=True,
                    levels=np.linspace(0,1,11), show_boundaries=False, show_mesh=True)
        mesh.refinements.future.h_adapt_to(f, levels=(0.5, 0.75))
        mesh = mesh.do.evolve()

    make_a_video_from_images_in_folder(image_folder, duration=8, clean_images=True)