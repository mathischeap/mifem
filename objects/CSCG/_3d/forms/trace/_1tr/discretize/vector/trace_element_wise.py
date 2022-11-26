# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('/')


from components.freeze.base import FrozenOnly
import numpy as np
from components.quadrature import Quadrature



class _3dCSCG_1Trace_Discretize_TraceElementWiseVector(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self.___cache_DISCRETIZE_STANDARD___ = None
        self.___cache_DISCRETIZE_TEW___ = None
        self._freeze_self_()

    def __call__(self,
        update_cochain=True, target='func', quad_degree=None):
        """We will discretize the Trace_parallel component of a standard vector field to all trace
        elements.


        'locally full local TEW cochain' means the cochain is a dict whose keys are trace-element
        numbers and values are trace-element-wise local cochains.

        """
        SELF = self._tf_
        # first check `target` and `update_cochain` inputs----------------------------------------
        if target in ('BC',):
            assert update_cochain is False, f"CANNOT update cochain when target is {target}"

        # ----- 3D: prepare and cache or read from cache the data for the numerical integration ----------------
        if self.___cache_DISCRETIZE_STANDARD___ is None or self.___cache_DISCRETIZE_STANDARD___['quadDegree'] != quad_degree:
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

            # NS ------------------------------------------------------
            qn_NS_dy_y = []
            qn_NS_dy_z = []
            for k in range(SELF.p[2]+1):
                for j in range(SELF.p[1]):
                    qn_NS_dy_y.append(qnodes[1][j])
                    qn_NS_dy_z.append(nodes[2][k]*np.ones(p[1]+1))
            qn_NS_dy_y, qn_NS_dy_z = np.array(qn_NS_dy_y), np.array(qn_NS_dy_z)
            lens_NS_dy = np.tile(lens[1]*0.5, (SELF.p[2]+1))

            qn_NS_dz_y = []
            qn_NS_dz_z = []
            for k in range(SELF.p[2]):
                for j in range(SELF.p[1]+1):
                    qn_NS_dz_y.append(nodes[1][j]*np.ones(p[2]+1))
                    qn_NS_dz_z.append(qnodes[2][k])
            qn_NS_dz_y, qn_NS_dz_z = np.array(qn_NS_dz_y), np.array(qn_NS_dz_z)
            lens_NS_dz = np.repeat(lens[2]*0.5, (SELF.p[1] + 1))

            # WE --------------------------------------------------------------------------------
            qn_WE_dx_x = []
            qn_WE_dx_z = []
            for k in range(SELF.p[2]+1):
                for i in range(SELF.p[0]):
                    qn_WE_dx_x.append(qnodes[0][i])
                    qn_WE_dx_z.append(nodes[2][k]*np.ones(p[0]+1))
            qn_WE_dx_x, qn_WE_dx_z = np.array(qn_WE_dx_x), np.array(qn_WE_dx_z)
            lens_WE_dx = np.tile(lens[0]*0.5, (SELF.p[2] + 1))


            qn_WE_dz_x = []
            qn_WE_dz_z = []
            for k in range(SELF.p[2]):
                for i in range(SELF.p[0]+1):
                    qn_WE_dz_x.append(nodes[0][i]*np.ones(p[2]+1))
                    qn_WE_dz_z.append(qnodes[2][k])
            qn_WE_dz_x, qn_WE_dz_z = np.array(qn_WE_dz_x), np.array(qn_WE_dz_z)
            lens_WE_dz = np.repeat(lens[2]*0.5, (SELF.p[0] + 1))

            #BF --------------------------------------------------------------------------------------------
            qn_BF_dx_x = []
            qn_BF_dx_y = []
            for j in range(SELF.p[1]+1):
                for i in range(SELF.p[0]):
                    qn_BF_dx_x.append(qnodes[0][i])
                    qn_BF_dx_y.append(nodes[1][j]*np.ones(p[0]+1))
            qn_BF_dx_x, qn_BF_dx_y = np.array(qn_BF_dx_x), np.array(qn_BF_dx_y)
            lens_BF_dx = np.tile(lens[0]*0.5, (SELF.p[1] + 1))

            qn_BF_dy_x = []
            qn_BF_dy_y = []
            for j in range(SELF.p[1]):
                for i in range(SELF.p[0]+1):
                    qn_BF_dy_x.append(nodes[0][i]*np.ones(p[1]+1))
                    qn_BF_dy_y.append(qnodes[1][j])
            qn_BF_dy_x, qn_BF_dy_y = np.array(qn_BF_dy_x), np.array(qn_BF_dy_y)
            lens_BF_dy = np.repeat(lens[1]*0.5, (SELF.p[0] + 1))

            LENS = [lens_NS_dy, lens_NS_dz, lens_WE_dx, lens_WE_dz, lens_BF_dx, lens_BF_dy]
            cd = dict()
            cd['quadDegree'] = quad_degree

            cd['LENS'] = LENS
            cd['qn_NS_dy_y'] = qn_NS_dy_y
            cd['qn_NS_dy_z'] = qn_NS_dy_z
            cd['qn_NS_dz_y'] = qn_NS_dz_y
            cd['qn_NS_dz_z'] = qn_NS_dz_z
            cd['qn_WE_dx_x'] = qn_WE_dx_x
            cd['qn_WE_dx_z'] = qn_WE_dx_z
            cd['qn_WE_dz_x'] = qn_WE_dz_x
            cd['qn_WE_dz_z'] = qn_WE_dz_z
            cd['qn_BF_dx_x'] = qn_BF_dx_x
            cd['qn_BF_dx_y'] = qn_BF_dx_y
            cd['qn_BF_dy_x'] = qn_BF_dy_x
            cd['qn_BF_dy_y'] = qn_BF_dy_y
            cd['quad_weights'] = quad_weights

            self.___cache_DISCRETIZE_STANDARD___ = cd
        else:
            qn_NS_dy_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_dy_y']
            qn_NS_dy_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_dy_z']
            qn_NS_dz_y = self.___cache_DISCRETIZE_STANDARD___['qn_NS_dz_y']
            qn_NS_dz_z = self.___cache_DISCRETIZE_STANDARD___['qn_NS_dz_z']
            qn_WE_dx_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_dx_x']
            qn_WE_dx_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_dx_z']
            qn_WE_dz_x = self.___cache_DISCRETIZE_STANDARD___['qn_WE_dz_x']
            qn_WE_dz_z = self.___cache_DISCRETIZE_STANDARD___['qn_WE_dz_z']
            qn_BF_dx_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_dx_x']
            qn_BF_dx_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_dx_y']
            qn_BF_dy_x = self.___cache_DISCRETIZE_STANDARD___['qn_BF_dy_x']
            qn_BF_dy_y = self.___cache_DISCRETIZE_STANDARD___['qn_BF_dy_y']

        # ----- 1D: prepare and cache or read from cache the data for the numerical integration ----------------
        if self.___cache_DISCRETIZE_TEW___ is None or self.___cache_DISCRETIZE_TEW___['quadDegree'] != quad_degree:
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
                qnodes.append(list(qnodes_i))

            # NS ------------------------------------------------------------------------------------
            qn_NS_dy_y1 = []
            qn_NS_dy_z1 = []
            for k in range(SELF.p[2]+1):
                for j in range(SELF.p[1]):
                    qn_NS_dy_y1.append(qnodes[1][j])
                    qn_NS_dy_z1.append([nodes[2][k],])
            lens_NS_dy = np.tile(lens[1]*0.5, (SELF.p[2]+1))

            qn_NS_dz_y1 = []
            qn_NS_dz_z1 = []
            for k in range(SELF.p[2]):
                for j in range(SELF.p[1]+1):
                    qn_NS_dz_y1.append([nodes[1][j],])
                    qn_NS_dz_z1.append(qnodes[2][k])
            lens_NS_dz = np.repeat(lens[2]*0.5, (SELF.p[1] + 1))

            # WE --------------------------------------------------------------------------------
            qn_WE_dx_x1 = []
            qn_WE_dx_z1 = []
            for k in range(SELF.p[2]+1):
                for i in range(SELF.p[0]):
                    qn_WE_dx_x1.append(qnodes[0][i])
                    qn_WE_dx_z1.append([nodes[2][k],])
            lens_WE_dx = np.tile(lens[0]*0.5, (SELF.p[2] + 1))


            qn_WE_dz_x1 = []
            qn_WE_dz_z1 = []
            for k in range(SELF.p[2]):
                for i in range(SELF.p[0]+1):
                    qn_WE_dz_x1.append([nodes[0][i],])
                    qn_WE_dz_z1.append(qnodes[2][k])
            lens_WE_dz = np.repeat(lens[2]*0.5, (SELF.p[0] + 1))

            #BF --------------------------------------------------------------------------------------------
            qn_BF_dx_x1 = []
            qn_BF_dx_y1 = []
            for j in range(SELF.p[1]+1):
                for i in range(SELF.p[0]):
                    qn_BF_dx_x1.append(qnodes[0][i])
                    qn_BF_dx_y1.append([nodes[1][j],])
            lens_BF_dx = np.tile(lens[0]*0.5, (SELF.p[1] + 1))

            qn_BF_dy_x1 = []
            qn_BF_dy_y1 = []
            for j in range(SELF.p[1]):
                for i in range(SELF.p[0]+1):
                    qn_BF_dy_x1.append([nodes[0][i],])
                    qn_BF_dy_y1.append(qnodes[1][j])
            lens_BF_dy = np.repeat(lens[1]*0.5, (SELF.p[0] + 1))

            LENS = [lens_NS_dy, lens_NS_dz, lens_WE_dx, lens_WE_dz, lens_BF_dx, lens_BF_dy]
            cd = dict()
            cd['quadDegree'] = quad_degree

            cd['LENS'] = LENS
            cd['qn_NS_dy_y1'] = qn_NS_dy_y1
            cd['qn_NS_dy_z1'] = qn_NS_dy_z1
            cd['qn_NS_dz_y1'] = qn_NS_dz_y1
            cd['qn_NS_dz_z1'] = qn_NS_dz_z1
            cd['qn_WE_dx_x1'] = qn_WE_dx_x1
            cd['qn_WE_dx_z1'] = qn_WE_dx_z1
            cd['qn_WE_dz_x1'] = qn_WE_dz_x1
            cd['qn_WE_dz_z1'] = qn_WE_dz_z1
            cd['qn_BF_dx_x1'] = qn_BF_dx_x1
            cd['qn_BF_dx_y1'] = qn_BF_dx_y1
            cd['qn_BF_dy_x1'] = qn_BF_dy_x1
            cd['qn_BF_dy_y1'] = qn_BF_dy_y1
            cd['quad_weights'] = quad_weights
            self.___cache_DISCRETIZE_TEW___ = cd
        else:
            LENS = self.___cache_DISCRETIZE_TEW___['LENS']
            qn_NS_dy_y1 = self.___cache_DISCRETIZE_TEW___['qn_NS_dy_y1']
            qn_NS_dy_z1 = self.___cache_DISCRETIZE_TEW___['qn_NS_dy_z1']
            qn_NS_dz_y1 = self.___cache_DISCRETIZE_TEW___['qn_NS_dz_y1']
            qn_NS_dz_z1 = self.___cache_DISCRETIZE_TEW___['qn_NS_dz_z1']
            qn_WE_dx_x1 = self.___cache_DISCRETIZE_TEW___['qn_WE_dx_x1']
            qn_WE_dx_z1 = self.___cache_DISCRETIZE_TEW___['qn_WE_dx_z1']
            qn_WE_dz_x1 = self.___cache_DISCRETIZE_TEW___['qn_WE_dz_x1']
            qn_WE_dz_z1 = self.___cache_DISCRETIZE_TEW___['qn_WE_dz_z1']
            qn_BF_dx_x1 = self.___cache_DISCRETIZE_TEW___['qn_BF_dx_x1']
            qn_BF_dx_y1 = self.___cache_DISCRETIZE_TEW___['qn_BF_dx_y1']
            qn_BF_dy_x1 = self.___cache_DISCRETIZE_TEW___['qn_BF_dy_x1']
            qn_BF_dy_y1 = self.___cache_DISCRETIZE_TEW___['qn_BF_dy_y1']
            quad_weights = self.___cache_DISCRETIZE_TEW___['quad_weights']

        #------- check func and get func --------------------------------------------------------------------
        if target == 'func':
            assert SELF.CF is not None, f"No func.body!"
            TEW_func = SELF.CF.___DO_evaluate_func_at_time___()
        elif target == 'BC':
            assert SELF.BC.CF is not None, f"No func.body!"
            TEW_func = SELF.BC.CF.___DO_evaluate_func_at_time___()
        else:
            raise NotImplementedError(f"1Trace = discretize_VectorField_standard of ftype trace-element-wise "
                                      f"not applicable for target={target}.")

        # --- Now, we do a check whether trace-elements in TEW_func are local -------------------------------
        for T in TEW_func:
            assert T in SELF.mesh.trace.elements, f"trace-element #{T} is not a local trace-element."

        # dispatch the lens for integration ------------------------------------------------
        lens_NS_dy, lens_NS_dz, lens_WE_dx, lens_WE_dz, lens_BF_dx, lens_BF_dy = LENS
        # initialize the Trace-Element-Wise-local-cochain dict ------------------------------------
        local_TEW = dict()
        Zo = [0,]
        # Go through all valid local trace-elements in the function ---------------------------
        for T in TEW_func:
            te = SELF.mesh.trace.elements[T]
            # ele = te.CHARACTERISTIC_element
            ele_side = te.CHARACTERISTIC_side

            #--------------------------------- NS sides ---------------------------------------
            if ele_side in 'NS':
                J = te.coordinate_transformation.Jacobian_matrix(qn_NS_dy_y, qn_NS_dy_z)
                J = (J[0][0], J[1][0], J[2][0]) # dy of (dy, dz)
                u, v, w = list(), list(), list()
                for dy, dz in zip(qn_NS_dy_y1, qn_NS_dy_z1):
                    ___, uvw = TEW_func[T](Zo, dy, dz)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,:,0]
                v = np.array(v)[:,:,0]
                w = np.array(w)[:,:,0]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[1]
                C = lens_NS_dy
                cochain_dy = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                J = te.coordinate_transformation.Jacobian_matrix(qn_NS_dz_y, qn_NS_dz_z)
                J = (J[0][1], J[1][1], J[2][1]) # dz of (dy, dz)
                u, v, w = list(), list(), list()
                for dy, dz in zip(qn_NS_dz_y1, qn_NS_dz_z1):
                    ___, uvw = TEW_func[T](Zo, dy, dz)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,0,:]
                v = np.array(v)[:,0,:]
                w = np.array(w)[:,0,:]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[2]
                C = lens_NS_dz
                cochain_dz = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dy, cochain_dz])

            #--------------------------------- WE sides ---------------------------------------
            elif ele_side in 'WE':
                J = te.coordinate_transformation.Jacobian_matrix(qn_WE_dx_x, qn_WE_dx_z)
                J = (J[0][1], J[1][1], J[2][1]) # dx of (dz, dx)
                u, v, w = list(), list(), list()
                for dx, dz in zip(qn_WE_dx_x1, qn_WE_dx_z1):
                    ___, uvw = TEW_func[T](dx, Zo, dz)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,:,0]
                v = np.array(v)[:,:,0]
                w = np.array(w)[:,:,0]

                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[0]
                C = lens_WE_dx
                cochain_dx = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                J = te.coordinate_transformation.Jacobian_matrix(qn_WE_dz_x, qn_WE_dz_z)
                J = (J[0][0], J[1][0], J[2][0]) # dz of (dz, dx)
                u, v, w = list(), list(), list()
                for dx, dz in zip(qn_WE_dz_x1, qn_WE_dz_z1):
                    ___, uvw = TEW_func[T](dx, Zo, dz)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,0,:]
                v = np.array(v)[:,0,:]
                w = np.array(w)[:,0,:]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[2]
                C = lens_WE_dz
                cochain_dz = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dx, cochain_dz])

            #--------------------------------- BF sides ---------------------------------------
            elif ele_side in 'BF':

                J = te.coordinate_transformation.Jacobian_matrix(qn_BF_dx_x, qn_BF_dx_y)
                J = (J[0][0], J[1][0], J[2][0]) # dx of (dx, dy)
                u, v, w = list(), list(), list()
                for dx, dy in zip(qn_BF_dx_x1, qn_BF_dx_y1):
                    ___, uvw = TEW_func[T](dx, dy, Zo)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,:,0]
                v = np.array(v)[:,:,0]
                w = np.array(w)[:,:,0]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[0]
                C = lens_BF_dx
                cochain_dx = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                J = te.coordinate_transformation.Jacobian_matrix(qn_BF_dy_x, qn_BF_dy_y)
                J = (J[0][1], J[1][1], J[2][1]) # dy of (dx, dy)
                u, v, w = list(), list(), list()
                for dx, dy in zip(qn_BF_dy_x1, qn_BF_dy_y1):
                    ___, uvw = TEW_func[T](dx, dy, Zo)
                    u.append(uvw[0])
                    v.append(uvw[1])
                    w.append(uvw[2])
                u = np.array(u)[:,0,:]
                v = np.array(v)[:,0,:]
                w = np.array(w)[:,0,:]
                A = J[0] * u + J[1] * v + J[2] * w
                B = quad_weights[1]
                C = lens_BF_dy
                cochain_dy = np.einsum('jk, k, j -> j', A, B, C, optimize='greedy')

                te_primal_local = np.concatenate([cochain_dx, cochain_dy])

            else:
                raise Exception()

            if not SELF.space.IS_Kronecker: raise NotImplementedError()

            local_TEW[T] = te_primal_local

        if update_cochain: SELF.cochain.local_TEW = local_TEW
        # 'locally full local TEW cochain': provide cochain.local_TEW and for all dofs on the trace element.
        return 'locally full local TEW cochain', local_TEW





if __name__ == '__main__':
    # mpiexec -n 5 python _3dCSCG\forms\trace\_1_trace\discretize\vector\trace_element_wise.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([2,2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)