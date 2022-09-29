# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/26 2:33 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

from root.config.main import rAnk, mAster_rank, cOmm, np
from screws.functions._3d_space.angle import angle_between_two_vectors

class _3dCSCG_MeshBoundaryCT_constant(FrozenOnly):
    """"""

    def __init__(self, boundary):
        """"""
        self._boundary_ = boundary
        self._constant_unit_normal_vector_ = True
        self._freeze_self_()

    @property
    def outward_unit_normal_vector(self):
        """If this boundary has a constant outward-unit-normal-vector, here we return it (a tuple of three
        number), otherwise, we return None.

        """
        if self._constant_unit_normal_vector_ is True: # DO NOT USE NONE FOR THIS ONE.

            boundary = self._boundary_

            local_trace_elements = boundary.trace_elements

            cV = None
            TF = False

            for i in local_trace_elements:
                trace_element = boundary.mesh.trace.elements[i]

                te_CV = trace_element.coordinate_transformation.constant.unit_normal_vector

                if te_CV is None:
                    TF = False
                    break
                else:
                    if cV is None:
                        cV = te_CV
                        TF = True
                    elif cV == te_CV:
                        TF = True
                    else:
                        TF = False
                        break
            if TF:
                LOCAL_CV = cV
            else:
                LOCAL_CV = None

            LOCAL_CV = cOmm.gather(LOCAL_CV, root=mAster_rank)

            if rAnk == mAster_rank:

                cN = LOCAL_CV.count(None)

                if cN == len(LOCAL_CV):
                    LOCAL_CV = None

                else:
                    for cv in LOCAL_CV:
                        if cv is not None:
                            break
                    LOCAL_CV = cv

            else:
                pass

            Constant_unit_normal_vector_ = cOmm.bcast(LOCAL_CV, root=mAster_rank)

            # now we need to check if this is outward, because for a trace element we don't know this!

            if Constant_unit_normal_vector_ is None:
                pass

            else:
                # this check may be wrong for some super strange shaped mesh elements.
                element_sides = boundary.element_sides

                _00_ = np.array([0, 0])
                CheckDict = {"N": [np.array([-0.75, -1]), _00_, _00_],
                             "S": [np.array([0.75, 1]), _00_, _00_],
                             "W": [_00_, np.array([-0.75, -1]), _00_],
                             "E": [_00_, np.array([0.75, 1]), _00_],
                             "B": [_00_, _00_, np.array([-0.75, -1])],
                             "F": [_00_, _00_, np.array([0.75, 1])],}
                As = list()
                for e_sd in element_sides:
                    e = int(e_sd[:-1])
                    sd = e_sd[-1]

                    xyz = CheckDict[sd]

                    element = boundary.mesh.elements[e]

                    eCT =element.coordinate_transformation
                    x, y, z = eCT.mapping(*xyz)
                    x = x[1] - x[0]
                    y = y[1] - y[0]
                    z = z[1] - z[0]

                    As.append(angle_between_two_vectors((x, y, z), Constant_unit_normal_vector_))

                if len(As) > 0:
                    oppo = all([_ > np.pi/2 for _ in As])
                    same = all([_ < np.pi/2 for _ in As])
                else:
                    oppo = None
                    same = None

                oppo = cOmm.gather(oppo, root=mAster_rank)
                same = cOmm.gather(same, root=mAster_rank)

                if rAnk == mAster_rank:

                    oppo = all([_ is True or _ is None for _ in oppo])
                    same = all([_ is True or _ is None for _ in same])

                    if oppo is False and same is False:
                        change = None
                    elif oppo:
                        assert not same
                        change = True
                    else:
                        assert same
                        change = False
                else:
                    change = None

                change = cOmm.bcast(change, root=mAster_rank)

                if change is None:
                    Constant_unit_normal_vector_ = None
                elif change:
                    a, b, c = Constant_unit_normal_vector_
                    Constant_unit_normal_vector_= (-a, -b, -c)
                else:
                    pass
            #============================================================================================

            self._constant_unit_normal_vector_ = Constant_unit_normal_vector_

        return self._constant_unit_normal_vector_

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
