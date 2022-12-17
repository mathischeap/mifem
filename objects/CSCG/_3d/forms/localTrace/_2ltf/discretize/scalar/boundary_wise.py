# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/1/2022 10:06 PM
"""
from components.freeze.main import FrozenOnly
import numpy as np
from components.quadrature import Quadrature


class _3dCSCG_2LocalTrace_DiscretizeBoundaryWise(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self.___cache_DISCRETIZE_STANDARD___ = None
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='func'):
        """

        Parameters
        ----------
        update_cochain
        target

        Returns
        -------

        """
        if target == 'func':
            FUNC = self._ltf_.CF.do.evaluate_func_at_time()
        elif target == 'BC':
            FUNC = self._ltf_.BC.CF.do.evaluate_func_at_time()
            assert update_cochain is False, \
                f"When target is {target}, cannot update cochain!"
        else:
            raise Exception()

        if self.___cache_DISCRETIZE_STANDARD___ is None:
            p = [self._ltf_.dqp[i] + 1 for i in range(self._ltf_.ndim)]

            quad_nodes, quad_weights = Quadrature(p, category='Gauss').quad
            nodes = self._ltf_.space.nodes
            num_edges = [len(nodes[i])-1 for i in range(self._ltf_.ndim)]
            lens = [nodes[i][1:]-nodes[i][0:-1] for i in range(self._ltf_.ndim)]
            qnodes = []
            for i in range(self._ltf_.ndim):
                qnodes_i = ((np.array(quad_nodes[i])+1)/2)[np.newaxis, :].repeat(
                    num_edges[i],
                    axis=0
                )*lens[i][:, np.newaxis]
                qnodes_i += np.array(nodes[i][:-1])[:, np.newaxis].repeat(p[i]+1, axis=1)
                qnodes.append(qnodes_i)
            # NS sides
            qn_NS_y = []
            qn_NS_z = []
            for k in range(self._ltf_.p[2]):
                for j in range(self._ltf_.p[1]):
                    qn_NS_y.append(qnodes[1][j][:, np.newaxis].repeat(p[2]+1, axis=1))
                    qn_NS_z.append(qnodes[2][k][np.newaxis, :].repeat(p[1]+1, axis=0))
            qn_NS_y, qn_NS_z = np.array(qn_NS_y), np.array(qn_NS_z)
            area_NS = np.kron(lens[2], lens[1]) * 0.25
            # WE sides
            qn_WE_x = []
            qn_WE_z = []
            for k in range(self._ltf_.p[2]):
                for i in range(self._ltf_.p[0]):
                    qn_WE_x.append(qnodes[0][i][:, np.newaxis].repeat(p[2]+1, axis=1))
                    qn_WE_z.append(qnodes[2][k][np.newaxis, :].repeat(p[0]+1, axis=0))
            qn_WE_x, qn_WE_z = np.array(qn_WE_x), np.array(qn_WE_z)
            area_WE = np.kron(lens[2], lens[0]) * 0.25
            # BF sides
            qn_BF_x = []
            qn_BF_y = []
            for j in range(self._ltf_.p[1]):
                for i in range(self._ltf_.p[0]):
                    qn_BF_x.append(qnodes[0][i][:, np.newaxis].repeat(p[1]+1, axis=1))
                    qn_BF_y.append(qnodes[1][j][np.newaxis, :].repeat(p[0]+1, axis=0))
            qn_BF_x, qn_BF_y = np.array(qn_BF_x), np.array(qn_BF_y)
            area_BF = np.kron(lens[1], lens[0]) * 0.25

            cd = dict()
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

        RANGE_element_sides = self._ltf_.mesh.boundaries.range_of_element_sides
        cochainLocal = dict()
        Tmap = self._ltf_.mesh.trace.elements.map
        NSWEBF = 'NSWEBF'
        for bn in FUNC:
            func_bn = FUNC[bn][0]
            element_sides = RANGE_element_sides[bn]
            elements, sides = list(), list()
            for element_side in element_sides:
                element = int(element_side[:-1])
                side = element_side[-1]
                elements.append(element)
                sides.append(side)

            for i, side in zip(elements, sides):
                j = Tmap[i][NSWEBF.index(side)]
                te = self._ltf_.mesh.trace.elements[j]

                if side in 'NS':
                    qn0, qn1 = qn_NS_y, qn_NS_z
                    qw0, qw1 = quad_weights[1], quad_weights[2]
                    area = area_NS
                elif side in 'WE':
                    qn0, qn1 = qn_WE_x, qn_WE_z
                    qw0, qw1 = quad_weights[0], quad_weights[2]
                    area = area_WE
                elif side in 'BF':
                    qn0, qn1 = qn_BF_x, qn_BF_y
                    qw0, qw1 = quad_weights[0], quad_weights[1]
                    area = area_BF
                else:
                    raise Exception()

                x, y, z = te.coordinate_transformation.mapping(
                    qn0, qn1, from_element=i, side=side
                )
                f = func_bn(x, y, z)
                sqrt_g = np.sqrt(
                    te.coordinate_transformation.metric(qn0, qn1)
                )
                cochain_i_side = np.einsum(
                    'mij, i, j, m -> m',
                    f * sqrt_g,
                    qw0, qw1, area,
                    optimize='greedy'
                )

                if i not in cochainLocal:
                    cochainLocal[i] = dict()
                else:
                    pass

                cochainLocal[i][side] = cochain_i_side

        if update_cochain:
            self._ltf_.cochain.local_ESW = cochainLocal
        else:
            pass

        return 'mesh-element-side-wise local cochain', cochainLocal
