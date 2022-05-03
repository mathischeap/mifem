# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import rAnk, mAster_rank, cOmm, np
from screws.freeze.main import FrozenOnly
import matplotlib.pyplot as plt

from screws.functions._3d_space.angle import angle_between_two_vectors
from itertools import combinations

from objects.CSCG._3d.mesh.trace.elements.do.find import _3dCSCG_Trace_Elements_Do_Find

class _3dCSCG_Trace_Elements_DO(FrozenOnly):
    """We find some specific groups of elements."""

    def __init__(self, trace_elements):
        self._elements_ = trace_elements
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _3dCSCG_Trace_Elements_Do_Find(self._elements_)
        return self._find_


    def illustrate_element(self, i, density_factor=2):
        """We illustrate the trace element #i with matplotlib on the
        mesh elements which share this trace element.

        To call this method, call it from all cores, otherwise, it does
        not work.

        :param int i: The trace element #i will be illustrated.
        :param int density_factor: be in [1,10], not be too large.
        :return:
        """
        # wo first find which core is going to do the plot _____________

        SELF = self._elements_

        if i in SELF:
            in_or_out = rAnk
            I_am_in = True
        else:
            in_or_out = -1
            I_am_in = False

        who_are_in = cOmm.gather(in_or_out, root=mAster_rank)
        if rAnk != mAster_rank:
            who_will_do_it = None
            in_how_many_cores = None
            the_other_core = None
        else:
            who_will_do_it = max(who_are_in)
            in_how_many_cores = 0
            the_other_core = None
            for _ in who_are_in:
                if _ != -1:
                    if _ !=who_will_do_it:
                        assert the_other_core is None
                        the_other_core = _
                    in_how_many_cores += 1

        who_will_do_it = cOmm.bcast(who_will_do_it, root=mAster_rank)
        in_how_many_cores = cOmm.bcast(in_how_many_cores, root=mAster_rank)
        the_other_core = cOmm.bcast(the_other_core, root=mAster_rank)

        sent_to = recv_from = None

        if in_how_many_cores == 1:
            if rAnk == who_will_do_it:
                assert SELF[i].IS_shared_by_cores is False, "trivial check"
            else:
                assert I_am_in is False, "trivial check"
        elif in_how_many_cores == 2:
            if I_am_in and rAnk != who_will_do_it:
                assert rAnk == the_other_core, "trivial check"
                # print(rAnk, ': I am going to provide data to rAnk ', who_will_do_it)
                sent_to = who_will_do_it
            elif I_am_in and rAnk == who_will_do_it:
                recv_from = the_other_core
            else:
                I_am_in = False
        else:
            raise Exception(f"can only have two cores sharing one trace element!")

        DATA2 = None
        r = s = cs = cr = uv_r = uv_s = None
        # prepare coordinate data first _____________________________________________
        if rAnk == who_will_do_it or I_am_in:
            density = 5 + 4 * density_factor
            i0 = 1 + density_factor
            i1 = 2 * density_factor + 2
            i2 = 3 * density_factor + 3
            _ = np.linspace(-1, 1, density)
            r, s = np.meshgrid(_, _, indexing='ij')
            anchors = ( # the points we will plot the outward unit norm vector
                [i0, i0],
                [i0, i2],
                [i2, i0],
                [i2, i2],
                [i1, i1],
            )
            uv_r = np.array([r[indices[0],indices[1]] for indices in anchors])
            uv_s = np.array([s[indices[0],indices[1]] for indices in anchors])
            cr, cs = np.array([r[i1, i1],]), np.array([s[i1, i1],])
            np.testing.assert_almost_equal(cr, 0)
            np.testing.assert_almost_equal(cs, 0)

        # Now let's prepare the data for the second subplot______________
        if in_how_many_cores == 1:
            if rAnk == who_will_do_it:
                if SELF[i].IS_on_mesh_boundary:
                    DATA2 = SELF[i].NON_CHARACTERISTIC_position # a string
                else:
                    ncp = SELF[i].NON_CHARACTERISTIC_position
                    other_element = int(ncp[:-1])
                    element = other_element
                    side = ncp[-1]
                    tes = SELF.map[element]
                    DATA2 = dict()
                    DATA2['element'] = element
                    DATA2['side'] = side
                    DATA2['tes'] = tes
                    for _, tei in enumerate(tes):
                        _side_ = 'NSWEBF'[_]
                        te = SELF._elements_[tei]
                        X, Y, Z = te.coordinate_transformation.mapping(
                            r, s, from_element=element, side=_side_)

                        if tei == i and _side_ == side:
                            x, y, z = te.coordinate_transformation.mapping(
                                uv_r, uv_s, from_element=element, side=_side_)
                            uvx, uvy, uvz = te.coordinate_transformation.\
                                ___PRIVATE_outward_unit_normal_vector___(
                                uv_r, uv_s, from_element=element, side=_side_)
                            DATA2[str(tei)+_side_] = [(X, Y, Z), (x, y, z), (uvx, uvy, uvz)]
                        else:
                            x, y, z = te.coordinate_transformation.mapping(
                                cr, cs, from_element=element, side=_side_)
                            DATA2[str(tei)+_side_] = [(X, Y, Z), (x, y, z), _side_]

            else:
                pass
        else: # in_how_many_cores == 2
            if rAnk == the_other_core:
                element = SELF[i].CHARACTERISTIC_element
                side = SELF[i].CHARACTERISTIC_side
                tes = SELF.map[element]
                DATA2 = dict()
                DATA2['element'] = element
                DATA2['side'] = side
                DATA2['tes'] = tes
                # noinspection PyTypeChecker
                assert tes.count(i) == 1, f"must be the case."
                for _, tei in enumerate(tes):
                    _side_ = 'NSWEBF'[_]

                    te = SELF._elements_[tei]
                    X, Y, Z = te.coordinate_transformation.mapping(
                        r, s, from_element=element, side=_side_)
                    if tei == i:
                        x, y, z = te.coordinate_transformation.mapping(
                            uv_r, uv_s, from_element=element, side=_side_)
                        uvx, uvy, uvz = te.coordinate_transformation.\
                            ___PRIVATE_outward_unit_normal_vector___(
                            uv_r, uv_s, from_element=element, side=_side_)
                        DATA2[str(tei)+_side_] = [(X, Y, Z), (x, y, z), (uvx, uvy, uvz)]
                    else:
                        x, y, z = te.coordinate_transformation.mapping(
                            cr, cs, from_element=element, side=_side_)
                        DATA2[str(tei)+_side_] = [(X, Y, Z), (x, y, z),  _side_   ]

                DATA2['rank'] = rAnk
                cOmm.send(DATA2, dest=sent_to, tag=who_will_do_it)
            elif rAnk == who_will_do_it:
                DATA2 = cOmm.recv(source=recv_from, tag=who_will_do_it)
            else:
                pass

        # let's do the plots ____________________________________________
        if rAnk == who_will_do_it:
            element = SELF[i].CHARACTERISTIC_element
            side = SELF[i].CHARACTERISTIC_side
            tes = SELF.map[element]
            assert i in tes, "trivial check"
            # noinspection PyTypeChecker
            assert 'NSWEBF'[tes.index(i)] == side

            fig = plt.figure(figsize=(14,6))
            TITLE = f"**trace element {i}**"
            if SELF[i].IS_on_periodic_boundary:
                TITLE += ' (periodic)'
            fig.suptitle(TITLE)

            ax = fig.add_subplot(121, projection='3d')
            for _, tei in enumerate(tes):
                te = SELF._elements_[tei]
                _side_ = 'NSWEBF'[_]
                X, Y, Z = te.coordinate_transformation.mapping(r, s, from_element=element, side=_side_)
                if tei == i and _side_ == side:
                    ax.plot_surface(X, Y, Z) # plot the trace element
                    x, y, z = te.coordinate_transformation.mapping(
                        uv_r, uv_s, from_element=element, side=_side_)
                    uvx, uvy, uvz = te.coordinate_transformation.\
                        ___PRIVATE_outward_unit_normal_vector___(
                        uv_r, uv_s, from_element=element, side=_side_)

                    x_range = np.max(X) - np.min(X)
                    y_range = np.max(Y) - np.min(Y)
                    z_range = np.max(Z) - np.min(Z)
                    mean_range = (x_range + y_range + z_range) / 6

                    ax.quiver(x, y, z, uvx*mean_range, uvy*mean_range, uvz*mean_range,
                              color='r', linewidth=0.5)
                    ax.text(x[-1] + 0.5*uvx[-1]*mean_range, y[-1] +
                            0.5*uvy[-1]*mean_range, z[-1] + 0.5*uvz[-1]*mean_range,
                            side, color='green', ha='center', va='center', ma='center')
                else:
                    ax.plot_surface(X, Y, Z, color=(0.7,0.7,0.7,0.5))
                    x, y, z = te.coordinate_transformation.mapping(
                        cr, cs, from_element=element, side=_side_)
                    ax.text(x[0], y[0], z[0], 'NSWEBF'[_],
                            color='k', ha='center', va='center', ma='center')

            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')

            if in_how_many_cores == 1 and isinstance(DATA2, str):
                plt.title(f"in rank {rAnk}, on [{side}] of element {element}.")
            elif in_how_many_cores == 1 and isinstance(DATA2, dict):
                plt.title(f"in rank {rAnk}, on [{side}] of characteristic element {element}.")
            elif in_how_many_cores == 2:
                plt.title(f"in rank {rAnk}, on [{side}] of characteristic element {element}.")
            else:
                raise Exception()

            ax = fig.add_subplot(122, projection='3d')
            if isinstance(DATA2, dict):
                element = DATA2['element']
                side_0 = side
                side = DATA2['side']
                tes = DATA2['tes']
                assert i in tes, "trivial check"
                # noinspection PyTypeChecker
                # assert 'NSWEBF'[tes.index(i)] == side
                assert side+side_0 in ('NS', 'SN', 'WE', 'EW', 'BF', 'FB')
                for _, tei in enumerate(tes):
                    _side_ = 'NSWEBF'[_]
                    if tei == i and _side_ == side:
                        XYZ, xyz, uv_xyz = DATA2[str(tei)+_side_]
                        ax.plot_surface(*XYZ)  # plot the trace element
                        uvx, uvy, uvz = uv_xyz
                        x, y, z = xyz
                        X, Y, Z = XYZ
                        x_range = np.max(X) - np.min(X)
                        y_range = np.max(Y) - np.min(Y)
                        z_range = np.max(Z) - np.min(Z)
                        mean_range = (x_range + y_range + z_range) / 6
                        ax.quiver(*xyz, uvx*mean_range, uvy*mean_range, uvz*mean_range,
                                  color='r', linewidth=0.5)
                        ax.text(x[-1] + 0.5*uvx[-1]*mean_range, y[-1] +
                                0.5*uvy[-1]*mean_range, z[-1] + 0.5*uvz[-1]*mean_range,
                                side, color='green', ha='center', va='center', ma='center')
                    else:
                        XYZ, xyz, _s_ = DATA2[str(tei)+_side_]
                        ax.plot_surface(*XYZ, color=(0.7, 0.7, 0.7, 0.5))
                        x, y, z = xyz
                        ax.text(x[0], y[0], z[0], _s_,
                                color='k', ha='center', va='center', ma='center')

            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_zlabel(r'$z$')
            if in_how_many_cores == 1 and isinstance(DATA2, str):
                plt.title(f"on mesh boundary: <{DATA2}>")
            elif in_how_many_cores == 1 and isinstance(DATA2, dict):
                plt.title(f"in rank {rAnk}, on [{side}] of element {element}.")
            elif in_how_many_cores == 2:
                plt.title(f"in rank {DATA2['rank']}, on [{side}] of characteristic element {element}.")
            else:
                raise Exception()
            plt.show()
            plt.close()
        else:
            return



    def get_quality_of_trace_elements(self):
        """

        :return: A dict whose keys are the numbers of local trace
            elements and values are the qualities of the trace elements.
        """
        SELF = self._elements_
        r = np.array([-0.5, -0.5, 0.5, 0.5, 0, 0.75, 0.75, -0.75, -0.75])
        s = np.array([-0.5, 0.5, -0.5, 0.5, 0, -0.75, 0.75, -0.75, 0.75])

        xi, et, sg = np.array([0, ]), np.array([0, ]), np.array([0, ])
        rc, sc = np.array([0, ]), np.array([0, ])

        LEN = len(r)
        assert len(s) == LEN, f"r, s length dis-match."
        Q = dict()
        Qs = list()
        for i in SELF:
            te = SELF[i]
            e = te.CHARACTERISTIC_element
            side = te.CHARACTERISTIC_side

            uv = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(
                r, s, from_element=e, side=side)
            vx, vy, vz = uv
            V = [[vx[_], vy[_], vz[_]] for _ in range(LEN)]
            comb = combinations(V, 2)
            angle = list()
            for v1v2 in comb:
                angle.append(angle_between_two_vectors(*v1v2))
            angle = max(angle)
            quality_0 =  float(1 - angle / np.pi)

            me = SELF._mesh_.elements[e]
            me_ct = me.coordinate_transformation.mapping(xi, et, sg)
            te_ct = te.coordinate_transformation.mapping(rc, sc, from_element=e, side=side)
            te_uv = te.coordinate_transformation.___PRIVATE_outward_unit_normal_vector___(
                rc, sc, from_element=e, side=side)
            v_inner = (me_ct[0]-te_ct[0], me_ct[1]-te_ct[1], me_ct[2]-te_ct[2])
            v_outer = te_uv
            angle = angle_between_two_vectors(v_inner, v_outer)
            quality_1 = float(angle / np.pi)


            _ = min([quality_0, quality_1])
            Q[i] = _
            Qs.append(_)

        if len(Qs) > 0:
            AvQ = sum(Qs) / len(Qs)
            WorstQ = min(Qs)
            BestQ = max(Qs)
        else:
            AvQ = None
            WorstQ = None
            BestQ = None

        return Q, AvQ, WorstQ, BestQ



if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\mesh\trace\elements\do.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [3, 4, 2]
    mesh = MeshGenerator('crazy_periodic', c=0.3, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh.trace.elements.selfcheck.outward_unit_normal_vector()
    # Q = mesh.trace.elements.quality

    mesh.trace.elements.do.illustrate_element(1)