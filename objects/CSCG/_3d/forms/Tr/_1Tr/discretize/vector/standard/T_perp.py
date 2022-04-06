

from screws.freeze.base import FrozenOnly



class _3dCSCG_1Tr_Discretize_StandardVector_T_perp(FrozenOnly):
    """"""
    def __init__(self, Tr):
        self._Tr_ = Tr
        self._freeze_self_()

    def __call__(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the Trace_perpendicular component of a standard vector field to all trace
        elements.

        :param update_cochain:
        :param target:
        :param quad_degree:
        :return:
        """
        if target in ('BC',): assert update_cochain is False, \
            f"CANNOT update cochain when target is {target}"
        raise NotImplementedError(quad_degree)





if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\forms\trace\_1_trace\discretize\vector\standard\T_perp.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)