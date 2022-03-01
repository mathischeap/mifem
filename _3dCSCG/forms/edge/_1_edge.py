# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from abc import ABC



from _3dCSCG.forms.edge.base.main import _3dCSCG_Edge



class _3dCSCG_1Edge(_3dCSCG_Edge, ABC):
    """
    Edge 1-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-1-edge-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_edge_1form')
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()


    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()


    def ___TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 1edge FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1edge FUNC do not accept func {func_body.__class__}")

    def ___TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 1edge BC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1edge BC do not accept func {func_body.__class__}")


    def discretize(self, update_cochain=True, target='func'):
        """
        Discretize the current function (a scalar field) to cochain.

        It is actually a wrapper of multiple methods that discretize functions of different types (a scalar
        field can be defined and represented in different ways in `python`, right?).

        :param bool update_cochain: (`default`: ``True``) If we update cochain with the output? Sometimes we
            may do not want to do so since we just want to use this method do some external jobs.
        :param target:
        :return: The cochain.
        :rtype: Its type can be different according to the particular discretize method.
        """
        if target == 'func':
            if self.TW.func.body.__class__.__name__ == '_3dCSCG_ScalarField':
                if self.func.ftype == 'standard':
                    return self.___PRIVATE_discretize_standard_ftype___(update_cochain=update_cochain)
                else:
                    raise NotImplementedError(f"3dCSCG 1-edge cannot (target func) discretize _3dCSCG_ScalarField of ftype={self.func.ftype}")

            else:
                raise NotImplementedError(f'3dCSCG 1-edge can not (target func) discretize {self.TW.func.body.__class__}.')

        else:
            raise NotImplementedError(f"3dCSCG 1-edge cannot discretize while targeting at {target}.")




    def ___PRIVATE_discretize_standard_ftype___(self, update_cochain=True, target='func'):
        """Discretize the standard _3dCSCG_ScalarField to a 1-edge-form

        'locally full local EEW cochain' means the cochain is a dict whose keys are edge-element
        numbers and values are edge-element-wise local cochains.

        """















if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\edge\_1_edge.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    e1 = FC('1-e')

    print(e1.standard_properties.tags)