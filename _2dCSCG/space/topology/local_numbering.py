


from root.config import *
from SCREWS.frozen import FrozenOnly
import matplotlib.pyplot as plt
import matplotlib.patches as patches



class LocalNumbering(FrozenOnly):
    """"""
    def __init__(self, FS):
        """"""
        assert FS.ndim == 2, " <LocalNumbering> "
        self._FS_ = FS
        self._freeze_self_()



    @property
    def _0Form_Inner(self):
        p = [self._FS_.p[i]+1 for i in range(self._FS_.ndim)]
        _ln_ = (np.arange(self._FS_.num_basis._0Form_Inner[0]).reshape(*p, order='F'),)
        return _ln_

    @property
    def _1Form_Inner(self):
        _ln_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[j]+1 for j in range(self._FS_.ndim)]
            p[i] -= 1
            I = 0 if i ==0 else np.sum(self._FS_.num_basis._1Form_Inner[1][0:i])
            _ln_ += (np.arange(self._FS_.num_basis._1Form_Inner[1][i]).reshape(*p, order='F') + I,)
        return _ln_

    @property
    def _2Form_Inner(self):
        return (np.arange(self._FS_.num_basis._2Form_Inner[0]).reshape(*self._FS_.p, order='F'),)





    @property
    def _0Form_Outer(self):
        p = [self._FS_.p[i]+1 for i in range(self._FS_.ndim)]
        _ln_ = (np.arange(self._FS_.num_basis._0Form_Outer[0]).reshape(*p, order='F'),)
        return _ln_

    @property
    def _1Form_Outer(self):
        _ln_ = ()
        for i in range(self._FS_.ndim):
            p = [self._FS_.p[j] for j in range(self._FS_.ndim)]
            p[i] += 1
            I = 0 if i == 0 else np.sum(self._FS_.num_basis._1Form_Outer[1][0:i])
            _ln_ += (np.arange(self._FS_.num_basis._1Form_Outer[1][i]).reshape(*p, order='F') + I,)
        return _ln_

    @property
    def _2Form_Outer(self):
        return (np.arange(self._FS_.num_basis._2Form_Outer[0]).reshape(*self._FS_.p, order='F'),)





    @property
    def _0Trace_Inner(self):
        """ """
        p = self._FS_.p
        _local_ = {'U': (np.arange(p[1]+1),),
                   'D': (np.arange(p[1]+1),),
                   'L': (np.arange(p[0]+1),),
                   'R': (np.arange(p[0]+1),),}
        return _local_

    @property
    def _0Trace_Outer(self):
        """ """
        p = self._FS_.p
        _local_ = {'U': (np.arange(p[1]+1),),
                   'D': (np.arange(p[1]+1),),
                   'L': (np.arange(p[0]+1),),
                   'R': (np.arange(p[0]+1),),}
        return _local_

    @property
    def _1Trace_Inner(self):
        p = self._FS_.p
        _local_ = {'U': (np.arange(p[1]),),
                   'D': (np.arange(p[1]),),
                   'L': (np.arange(p[0]),),
                   'R': (np.arange(p[0]),)}
        return _local_

    @property
    def _1Trace_Outer(self):
        p = self._FS_.p
        _local_ = {'U': (np.arange(p[1]),),
                   'D': (np.arange(p[1]),),
                   'L': (np.arange(p[0]),),
                   'R': (np.arange(p[0]),)}
        return _local_




    def ___PRIVATE_matplot___(self, what, **kwargs):
        """"""
        if what[2:6] == 'Form': # form
            self.___PRIVATE_matplot_form___(what, **kwargs)
        else:
            raise NotImplementedError()

    def ___PRIVATE_matplot_form___(self, what_form, saveto=None, usetex=False):
        """"""
        if rAnk == mAster_rank:
            fig, ax = plt.subplots(figsize=(14, 10))
            if 'Inner' in what_form:
                LN0 = self._0Form_Inner
                LN1 = self._1Form_Inner
                LN2 = self._2Form_Inner
                orientation = 'Inner'
            elif 'Outer' in what_form:
                LN0 = self._0Form_Outer
                LN1 = self._1Form_Outer
                LN2 = self._2Form_Outer
                orientation = 'Outer'
            else:
                raise Exception()
            k = int(what_form[1])
            self.___PRIVATE_plot_reference_mesh___(ax, orientation, k, usetex=usetex)
            LN = [LN0, LN1, LN2]
            XI, ETA = self._FS_.nodes
            p = self._FS_.p
            for i in range(3):
                if i == k:
                    color = 'b'
                else:
                    color = 'k'
                ln = LN[i]
                if i == 0: # plot node dofs
                    for m in range(p[0]+1):
                        for n in range(p[1]+1):
                            x = XI[m]
                            y = ETA[n]
                            numbering = ln[0][m,n]
                            plt.text(x, y, numbering, c=color, va='center', ha='center')

                elif i == 1: # plot edges
                    if 'Inner' in what_form:
                        for m in range(p[0]):
                            for n in range(p[1]+1):
                                x = (XI[m+1] + XI[m])/2
                                y = ETA[n]
                                numbering = ln[0][m,n]
                                plt.text(x, y, numbering, c=color, va='center', ha='center')
                        for m in range(p[0]+1):
                            for n in range(p[1]):
                                x = XI[m]
                                y = (ETA[n+1] + ETA[n]) / 2
                                numbering = ln[1][m, n]
                                plt.text(x, y, numbering, c=color, va='center', ha='center')
                    elif 'Outer' in what_form:
                        for m in range(p[0]):
                            for n in range(p[1]+1):
                                x = (XI[m+1] + XI[m])/2
                                y = ETA[n]
                                numbering = ln[1][m,n]
                                plt.text(x, y, numbering, c=color, va='center', ha='center')
                        for m in range(p[0]+1):
                            for n in range(p[1]):
                                x = XI[m]
                                y = (ETA[n+1] + ETA[n]) / 2
                                numbering = ln[0][m, n]
                                plt.text(x, y, numbering, c=color, va='center', ha='center')
                    else:
                        raise Exception()

                elif i == 2: # plot faces
                    for m in range(p[0]):
                        for n in range(p[1]):
                            x = (XI[m+1] + XI[m])/2
                            y = (ETA[n+1] + ETA[n]) / 2
                            numbering = ln[0][m,n]
                            plt.text(x, y, numbering, c=color, va='center', ha='center')

                else:
                    raise Exception()

            plt.title(f'local numbering of {orientation} {k}-form in' + r' $\Omega_{\mathrm{ref}}$')

            plt.show()
            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')

            return fig

    def ___PRIVATE_plot_reference_mesh___(self, ax, orientation, k, usetex=False):
        """"""
        if rAnk == mAster_rank:
            plt.rc('text', usetex=usetex)
            ax.set_aspect('equal')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
            plt.xlabel(r'$\xi$')
            plt.ylabel(r'$\eta$')

            plt.plot([-1,1],[-1,-1], '--',c='k')
            plt.plot([-1,1],[1,1], '--', c='k')
            plt.plot([-1,-1],[-1,1], '--', c='k')
            plt.plot([1,1],[-1,1], '--', c='k')

            nodes_xi, nodes_eta = self._FS_.nodes
            for n in nodes_xi:
                plt.plot([n,n],[-1,1], c='k')
            for n in nodes_eta:
                plt.plot([-1,1],[n,n], c='k')
            p = self._FS_.p
            colors = ["lightseagreen", "lightseagreen", "lightseagreen"]
            colors[k] = 'red'

            # -----------------------------------------------------------------------------------
            if orientation == 'Inner':
                # inner nodes ...
                A = list()
                style = f"-|>,head_width=3,head_length=5"
                kw = dict(arrowstyle=style, color=colors[0])
                for i in range(p[0]+1):
                    for j in range(p[1]+1):

                        if i == 0:
                            edge_x = nodes_xi[1] - nodes_xi[0]
                        elif i == p[0]:
                            edge_x = nodes_xi[-1] - nodes_xi[-2]
                        else:
                            edge_x1 = nodes_xi[i+1] - nodes_xi[i]
                            edge_x2 = nodes_xi[i] - nodes_xi[i-1]
                            edge_x = min([edge_x1, edge_x2])

                        if j == 0:
                            edge_y = nodes_eta[1] - nodes_eta[0]
                        elif j == p[1]:
                            edge_y = nodes_eta[-1] - nodes_eta[-2]
                        else:
                            edge_y1 = nodes_eta[j+1] - nodes_eta[j]
                            edge_y2 = nodes_eta[j] - nodes_eta[j-1]
                            edge_y = min([edge_y1, edge_y2])

                        edge = min([edge_x, edge_y])

                        d = 0.18 * edge

                        p0, p1 = nodes_xi[i], nodes_eta[j]
                        A.append( patches.FancyArrowPatch((p0-d, p1-d), (p0, p1), **kw) )
                        A.append( patches.FancyArrowPatch((p0+d, p1+d), (p0, p1), **kw) )
                        A.append( patches.FancyArrowPatch((p0+d, p1-d), (p0, p1), **kw) )
                        A.append( patches.FancyArrowPatch((p0-d, p1+d), (p0, p1), **kw) )
                for a in A:
                    plt.gca().add_patch(a)

                # inner edges ...
                for i in range(p[0]+1):
                    for j in range(p[1]):
                        p0 = nodes_xi[i]
                        edge_y = nodes_eta[j+1] - nodes_eta[j]
                        p1 = (nodes_eta[j+1] + nodes_eta[j])/2
                        al = edge_y * 0.1
                        plt.plot([p0, p0-al/2], [p1+al/2, p1-al/2], color=colors[1])
                        plt.plot([p0, p0+al/2], [p1+al/2, p1-al/2], color=colors[1])
                for i in range(p[0]):
                    for j in range(p[1]+1):
                        p0 = (nodes_xi[i+1] + nodes_xi[i])/2
                        edge_x = nodes_xi[i+1] - nodes_xi[i]
                        p1 = nodes_eta[j]
                        al = edge_x * 0.1
                        plt.plot([p0+al/2, p0-al/2], [p1, p1+al/2], color=colors[1])
                        plt.plot([p0+al/2, p0-al/2], [p1, p1-al/2], color=colors[1])

                # inner faces ...
                A = list()
                for i in range(p[0]):
                    for j in range(p[1]):
                        edge_x = nodes_xi[i+1] - nodes_xi[i]
                        edge_y = nodes_eta[j+1] - nodes_eta[j]
                        edge = min([edge_x, edge_y])
                        rad = 0.6 * edge
                        center = [(nodes_xi[i+1] + nodes_xi[i])/2, (nodes_eta[j+1] + nodes_eta[j])/2]
                        # noinspection PyTypeChecker
                        A.append(patches.Arc(center, width=rad, height=rad,
                                             theta1=180, theta2=90, color=colors[2]))
                        al = rad/8
                        plt.plot([center[0], center[0]+al*0.75],
                                 [center[1]+rad/2, center[1]+rad/2-al*0.8], color=colors[2])
                        plt.plot([center[0], center[0]+al],
                                 [center[1]+rad/2, center[1]+rad/2+al/2], color=colors[2])
                for a in A:
                    plt.gca().add_patch(a)
            # -----------------------------------------------------------------------------------
            elif orientation == 'Outer':
                # outer nodes ...
                A = list()
                for i in range(p[0]+1):
                    for j in range(p[1]+1):
                        p0, p1 = nodes_xi[i], nodes_eta[j]

                        if i == 0:
                            edge_x = nodes_xi[1] - nodes_xi[0]
                        elif i == p[0]:
                            edge_x = nodes_xi[-1] - nodes_xi[-2]
                        else:
                            edge_x1 = nodes_xi[i+1] - nodes_xi[i]
                            edge_x2 = nodes_xi[i] - nodes_xi[i-1]
                            edge_x = min([edge_x1, edge_x2])

                        if j == 0:
                            edge_y = nodes_eta[1] - nodes_eta[0]
                        elif j == p[1]:
                            edge_y = nodes_eta[-1] - nodes_eta[-2]
                        else:
                            edge_y1 = nodes_eta[j+1] - nodes_eta[j]
                            edge_y2 = nodes_eta[j] - nodes_eta[j-1]
                            edge_y = min([edge_y1, edge_y2])

                        edge = min([edge_x, edge_y])
                        center = [p0, p1]
                        rad = 0.3 * edge
                        # noinspection PyTypeChecker
                        A.append(patches.Arc(center, width=rad, height=rad,
                                             theta1=180, theta2=90, color=colors[0]))
                        al = rad/8
                        plt.plot([center[0], center[0]+al*0.75],
                                 [center[1]+rad/2, center[1]+rad/2-al*0.8], color=colors[0])
                        plt.plot([center[0], center[0]+al],
                                 [center[1]+rad/2, center[1]+rad/2+al/2], color=colors[0])
                for a in A:
                    plt.gca().add_patch(a)

                # outer edges ...
                A = list()
                style = f"->,head_width=5,head_length=10"
                kw = dict(arrowstyle=style, color=colors[1])
                for i in range(p[0]+1):
                    for j in range(p[1]):
                        if i == 0:
                            edge_x = nodes_xi[1] - nodes_xi[0]
                        elif i == p[0]:
                            edge_x = nodes_xi[-1] - nodes_xi[-2]
                        else:
                            edge_x1 = nodes_xi[i+1] - nodes_xi[i]
                            edge_x2 = nodes_xi[i] - nodes_xi[i-1]
                            edge_x = min([edge_x1, edge_x2])
                        al = edge_x * 0.5
                        p0 = nodes_xi[i]
                        p1 = (nodes_eta[j+1] + nodes_eta[j])/2
                        A.append(patches.FancyArrowPatch((p0 - al/2, p1), (p0+al/2, p1), **kw))

                for i in range(p[0]):
                    for j in range(p[1]+1):
                        if j == 0:
                            edge_y = nodes_eta[1] - nodes_eta[0]
                        elif j == p[1]:
                            edge_y = nodes_eta[-1] - nodes_eta[-2]
                        else:
                            edge_y1 = nodes_eta[j+1] - nodes_eta[j]
                            edge_y2 = nodes_eta[j] - nodes_eta[j-1]
                            edge_y = min([edge_y1, edge_y2])
                        al = edge_y * 0.5
                        p0 = (nodes_xi[i+1] + nodes_xi[i])/2
                        p1 = nodes_eta[j]
                        A.append(patches.FancyArrowPatch((p0, p1-al/2), (p0, p1+al/2), **kw))

                for a in A:
                    plt.gca().add_patch(a)

                # outer faces ...
                A = list()
                style = f"<-,head_width=3,head_length=6"
                kw = dict(arrowstyle=style, color=colors[2])

                for i in range(p[0]):
                    for j in range(p[1]):
                        edge_x = nodes_xi[i+1] - nodes_xi[i]
                        edge_y = nodes_eta[j+1] - nodes_eta[j]
                        edge = min([edge_x, edge_y])

                        d = edge * 0.18

                        p0, p1 = (nodes_xi[i+1] + nodes_xi[i])/2, (nodes_eta[j+1] + nodes_eta[j])/2
                        A.append( patches.FancyArrowPatch((p0-d, p1-d), (p0, p1), **kw) )
                        A.append( patches.FancyArrowPatch((p0+d, p1+d), (p0, p1), **kw) )
                        A.append( patches.FancyArrowPatch((p0+d, p1-d), (p0, p1), **kw) )
                        A.append( patches.FancyArrowPatch((p0-d, p1+d), (p0, p1), **kw) )
                for a in A:
                    plt.gca().add_patch(a)