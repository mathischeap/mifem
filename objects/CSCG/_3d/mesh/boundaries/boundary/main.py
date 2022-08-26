"""


"""


import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.boundaries.boundary.visualize import _3dCSCG_Mesh_Boundary_VIS

from root.config.main import rAnk, mAster_rank, cOmm, np
from screws.functions._3d_space.angle import angle_between_two_vectors

class _3dCSCG_Mesh_Boundary(FrozenOnly):
    def __init__(self, bdrs, name):
        self._bdrs_ = bdrs
        self._name_ = name
        self._visualize_ = None
        self._constant_unit_normal_vector_ = True
        self._freeze_self_()

    @property
    def name(self):
        return self._name_

    @property
    def mesh(self):
        return self._bdrs_._mesh_

    @property
    def element_sides(self):
        """This mesh boundary covers these local mesh element sides."""
        return self._bdrs_.range_of_element_sides[self._name_]

    @property
    def trace_elements(self):
        """This mesh boundary covers these local trace elements."""
        return self._bdrs_.range_of_trace_elements[self._name_]

    @property
    def region_sides(self):
        """This mesh boundary globally covers these region sides."""
        return self._bdrs_.range_of_region_sides[self._name_]

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_Mesh_Boundary_VIS(self)
        return self._visualize_

    @property
    def constant_unit_normal_vector(self):
        """outward by default.

        If this boundary has a constant outward-unit-normal-vector, here we return it (a tuple of three
        number), otherwise, we return None.

        """
        if self._constant_unit_normal_vector_ is True: # we need to compute it, DO NOT USE NONE FOR THIS ONE.

            local_trace_elements = self.trace_elements

            cV = None
            TF = False

            for i in local_trace_elements:
                trace_element = self.mesh.trace.elements[i]

                te_CV = trace_element.constant_unit_normal_vector

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
                element_sides = self.element_sides

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

                    element = self.mesh.elements[e]

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


if __name__ == '__main__':
    # mpiexec -n 8 python _3dCSCG\mesh\boundaries\boundary\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [5, 5, 5]
    # mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    boundaries = mesh.boundaries
    boundary = boundaries['Bottom']

    boundary.visualize()