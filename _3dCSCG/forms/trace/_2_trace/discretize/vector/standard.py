
import sys
if './' not in sys.path: sys.path.append('./')


from screws.freeze.inheriting.frozen_only import FrozenOnly
import numpy as np
from screws.quadrature import Quadrature



class _3dCSCG_2Trace_Discretize_StandardVector(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self.___cache_DISCRETIZE_STANDARD___ = None
        self._freeze_self_()

    def __call__(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the norm component of the vector to the trace 2-form.

        As a trace 2-form can only accommodate scalar, when given a
        vector, for example denoted by vec(u), we actually mean we get
        vec(u) · vec(n) where vec(n) represent the positive normal
        direction of the surface (trace element). The positive normal
        direction means x+, y+ or z+ direction, or more accurately,
        xi+, eta+ or sigma+ direction (in the mesh element setting)
        before mapping.


        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.
        """
        SELF = self._tf_

        if target in ('BC',): assert update_cochain is False, f"CANNOT update cochain when target is {target}"

        if self.___cache_DISCRETIZE_STANDARD___ is None or \
            self.___cache_DISCRETIZE_STANDARD___['quadDegree'] != quad_degree:
            p = [SELF.dqp[i] + 1 for i in range(SELF.ndim)] if quad_degree is None else quad_degree
            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
            nodes = SELF.space.nodes
            num_edges = [len(nodes[i])-1 for i in range(SELF.ndim)]
            lens = [nodes[i][1:]-nodes[i][0:-1] for i in range(SELF.ndim)]
            qnodes = []
            for i in range(SELF.ndim):
                qnodes_i = ((np.array(quad_nodes[i])+1)/2)[np.newaxis,:].repeat(num_edges[i],
                           axis=0)*lens[i][:,np.newaxis]
                qnodes_i += np.array(nodes[i][:-1])[:,np.newaxis].repeat(p[i]+1, axis=1)
                qnodes.append(qnodes_i)
            # NS sides
            qn_NS_y = []
            qn_NS_z = []
            for k in range(SELF.p[2]):
                for j in range(SELF.p[1]):
                    qn_NS_y.append(qnodes[1][j][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_NS_z.append(qnodes[2][k][np.newaxis,:].repeat(p[1]+1, axis=0))
            qn_NS_y, qn_NS_z = np.array(qn_NS_y), np.array(qn_NS_z)
            area_NS = np.kron(lens[2], lens[1]) * 0.25
            # WE sides
            qn_WE_x = []
            qn_WE_z = []
            for k in range(SELF.p[2]):
                for i in range(SELF.p[0]):
                    qn_WE_x.append(qnodes[0][i][:,np.newaxis].repeat(p[2]+1, axis=1))
                    qn_WE_z.append(qnodes[2][k][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_WE_x, qn_WE_z = np.array(qn_WE_x), np.array(qn_WE_z)
            area_WE = np.kron(lens[2], lens[0]) * 0.25
            # BF sides
            qn_BF_x = []
            qn_BF_y = []
            for j in range(SELF.p[1]):
                for i in range(SELF.p[0]):
                    qn_BF_x.append(qnodes[0][i][:,np.newaxis].repeat(p[1]+1, axis=1))
                    qn_BF_y.append(qnodes[1][j][np.newaxis,:].repeat(p[0]+1, axis=0))
            qn_BF_x, qn_BF_y = np.array(qn_BF_x), np.array(qn_BF_y)
            area_BF = np.kron(lens[1], lens[0]) * 0.25

            cd = dict()
            cd['quadDegree'] = quad_degree
            cd['qn_NS_y'] = qn_NS_y
            cd['qn_NS_z'] = qn_NS_z
            cd['area_NS'] = area_NS
            cd['qn_WE_x'] = qn_WE_x
            cd['qn_WE_z'] = qn_WE_z
            cd['area_WE'] = area_WE
            cd['qn_BF_x'] = qn_BF_x
            cd['qn_BF_y'] = qn_BF_y
            cd['area_BF'] = area_BF
            cd['quad_weights'] = quad_weights
            self.___cache_DISCRETIZE_STANDARD___ = cd
        else:
            qn_NS_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_y']
            qn_NS_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_z']
            area_NS = self.___cache_DISCRETIZE_STANDARD___['area_NS']
            qn_WE_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_x']
            qn_WE_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_z']
            area_WE = self.___cache_DISCRETIZE_STANDARD___['area_WE']
            qn_BF_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_x']
            qn_BF_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_y']
            area_BF = self.___cache_DISCRETIZE_STANDARD___['area_BF']
            quad_weights = self.___cache_DISCRETIZE_STANDARD___['quad_weights']

        if target == 'func':
            assert SELF.func.body is not None, f"No func.body!"
            _lf0_, _lf1_, _lf2_ = SELF.func.body
        elif target == 'BC':
            assert SELF.BC.body is not None, f"No BC.body!"
            _lf0_, _lf1_, _lf2_ = SELF.BC.body
        else:
            raise NotImplementedError(f"Not applicable for target={target}.")

        local_TEW = dict()
        for key in SELF.mesh.trace.elements:
            te = SELF.mesh.trace.elements[key]
            ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side
            if ele_side in 'NS':
                qn0, qn1 = qn_NS_y, qn_NS_z
                qw0, qw1 = quad_weights[1], quad_weights[2]
                area = area_NS
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
                ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
            elif ele_side in 'WE':
                qn0, qn1 = qn_WE_x, qn_WE_z
                qw0, qw1 = quad_weights[0], quad_weights[2]
                area = area_WE
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
                ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
            elif ele_side in 'BF':
                qn0, qn1 = qn_BF_x, qn_BF_y
                qw0, qw1 = quad_weights[0], quad_weights[1]
                area = area_BF
                x, y, z = te.coordinate_transformation.mapping(qn0, qn1, from_element=ele, side=ele_side)
                g = te.coordinate_transformation.metric(qn0, qn1)
                ouv = te.coordinate_transformation.unit_normal_vector(qn0, qn1)
            else:
                raise Exception()

            f0, f1, f2 = _lf0_(x, y, z), _lf1_(x, y, z), _lf2_(x, y, z)
            f = f0 * ouv[0] + f1 * ouv[1] + f2 * ouv[2]
            sqrt_g = np.sqrt(g)
            te_primal_local = np.einsum('mij, i, j, m -> m', f * sqrt_g,
                                        qw0, qw1, area,
                                        optimize='greedy')

            if not SELF.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[key] = te_primal_local


        if update_cochain: SELF.cochain.local_TEW = local_TEW
        # 'locally full local TEW cochain': provide cochain.local_TEW and for all dofs on the trace element.
        return 'locally full local TEW cochain', local_TEW





if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\forms\trace\_2_trace\discretize\scalar\vector.py

    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)