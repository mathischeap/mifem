# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')

from root.config.main import *
from components.freeze.main import FrozenOnly
import matplotlib.pyplot as plt




class _3dCSCG_standard_form_Matplot(FrozenOnly):
    """"""

    def __init__(self, sf):
        """ """
        assert '3dCSCG_standard_form' in sf.standard_properties.tags
        self._sf_ = sf

    def perpendicular_slice(self, MPS, plot_type='contourf', usetex=True, colormap='coolwarm',
        numOfSamples=100000, figsize=(6, 5),
        num_of_levels=20,
        xlabel=None, ylabel=None,
        title=None, levels=None, # if provide them, put them in list of length 1 (for 0-, 3-form) or 3 (for 1-, 2-form)
        colorbar_font_size=12, title_pad=10,
        label_size=12, title_size=12,
        minor_tick_length=4, major_tick_length=7, tick_pad=8, tick_size=12,
        saveto=None,):
        """

        :param MPS: A _3dCSCG_MeshPerpendicularSlice instance.
        :param plot_type:
        :param figsize:
        :param numOfSamples:
        :param usetex:
        :param colormap:
        :param num_of_levels: Only valid when ``levels`` is ``None``.
        :param xlabel:
        :param ylabel:
        :param title:
        :param levels:
        :param colorbar_font_size:
        :param title_pad:
        :param label_size:
        :param title_size:
        :param minor_tick_length:
        :param major_tick_length:
        :param tick_pad:
        :param tick_size:
        :param saveto:
        :return:
        """
        assert MPS.__class__.__name__ == "_3dCSCG_MeshPerpendicularSlice", "I need a _3dCSCG_MeshPerpendicularSlice."
        assert MPS._mesh_ == self._sf_.mesh, "Meshes do not match."
        assert self._sf_.cochain._local_ is not None, "Form has no cochain!"

        PTA = MPS.perpendicular_to_axis

        loc_len = len(MPS)
        ALL_LEN = COMM.gather(loc_len, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ALL_LEN = sum(ALL_LEN)
            assert ALL_LEN > 0, "MPS is empty!"
            density = int(np.sqrt(numOfSamples / ALL_LEN)) + 1
        else:
            density = None
        density = COMM.bcast(density, root=MASTER_RANK)
        sample = np.linspace(-1, 1, density)

        XYZ, VAL = dict(), dict() # get data for plot.
        for e in MPS:
            eps = MPS[e]
            pta = eps.perpendicular_to_axis
            assert pta == PTA, "For _3dCSCG_MeshPerpendicularSlice, all element perpendicular slice must have same PTA."
            pos = eps.position
            if pta == 'xi':
                xi = np.array([pos,])
                eta = sigma = sample
            elif pta == 'eta':
                eta = np.array([pos,])
                xi = sigma = sample
            elif pta == 'sigma':
                sigma = np.array([pos,])
                xi = eta = sample
            else:
                raise Exception()

            xyz, val = self._sf_.reconstruct(xi, eta, sigma, ravel=False, i=eps._element_.i)
            val = val[e]
            xyz = xyz[e]

            if self._sf_.k in (0, 3):
                assert len(val) == 1
            else:
                assert len(val) == 3

            _val_ = list()
            _xyz_ = list()

            for i, xyz_i in enumerate(xyz):
                if pta == 'xi':
                    _xyz_.append(xyz_i[0,:,:])
                elif pta == 'eta':
                    _xyz_.append(xyz_i[:,0,:])
                elif pta == 'sigma':
                    _xyz_.append(xyz_i[:,:,0])
                else:
                    raise Exception()

            for i, vi in enumerate(val):
                if pta == 'xi':
                    _val_.append(vi[0,:,:])
                elif pta == 'eta':
                    _val_.append(vi[:,0,:])
                elif pta == 'sigma':
                    _val_.append(vi[:,:,0])
                else:
                    raise Exception()

            XYZ[e] = _xyz_
            VAL[e] = _val_

        # gather all information to the mAster core ----------- BELOW ---------------------------------------------
        XYZ = COMM.gather(XYZ, root=MASTER_RANK)
        VAL = COMM.gather(VAL, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            ___ = dict()
            for xyz in XYZ:
                ___.update(xyz)
            XYZ = ___

            ___ = dict()
            for val in VAL:
                ___.update(val)
            VAL = ___

            del ___


            if self._sf_.k in (0, 3):
                MIN, MAX = [None,], [None,]

            else:
                MIN, MAX = [None, None, None], [None, None, None]

            for e in VAL:
                for i, val in enumerate(VAL[e]):
                    min_i = np.min(val)
                    max_i = np.max(val)

                    if MIN[i] is None:
                        MIN[i] = min_i
                    else:
                        if min_i < MIN[i]:
                            MIN[i] = min_i
                        else:
                            pass

                    if MAX[i] is None:
                        MAX[i] = max_i
                    else:
                        if max_i > MAX[i]:
                            MAX[i] = max_i
                        else:
                            pass

            for i, _ in enumerate(MIN):
                if MIN[i] == MAX[i]:
                    MIN[i] -= 0.5
                    MAX[i] += 0.5

            if levels is None:
                levels = list()
                for i, _ in enumerate(MIN):
                    levels.append(np.linspace(MIN[i], MAX[i], num_of_levels))


        # Now, we can do the plot ------------- BELOW -----------------------------------------------------------

        if RANK == MASTER_RANK:

            if self._sf_.k in (0, 3):
                NUM_PLOT = 1
                if title is None:
                    title = [f'{self._sf_.k}-form: {self._sf_.standard_properties.name}',]
            else:
                NUM_PLOT = 3
                if title is None:
                    title = list()
                    for _ in range(3):
                        title.append(f'{_}th component of {self._sf_.k}-form: {self._sf_.standard_properties.name}')

            plt.rc('text', usetex=usetex)
            if colormap is not None: plt.rcParams['image.cmap'] = colormap

            plotter = getattr(plt, plot_type)
            FIGURES = list()
            for n in range(NUM_PLOT):

                fig = plt.figure(figsize=figsize)

                for e in VAL: # go through all involved elements.
                    if PTA == 'xi':
                        axis_1, axis_2 = XYZ[e][1], XYZ[e][2]
                    elif PTA == 'eta':
                        axis_1, axis_2 = XYZ[e][0], XYZ[e][2]
                    elif PTA == 'sigma':
                        axis_1, axis_2 = XYZ[e][0], XYZ[e][1]
                    else:
                        raise Exception()

                    plotter(axis_1, axis_2, VAL[e][n], levels=levels[n])

                cb = plt.colorbar()
                # cb.set_label("Relative Photon Intensity", labelpad=-1, size=14) # change color bar name, gap and size.
                cb.ax.tick_params(labelsize=colorbar_font_size)

                if xlabel is None:
                    if PTA == 'xi':
                        xlabel = r"$\eta$"
                    elif PTA == 'eta':
                        xlabel = r"$\xi$"
                    elif PTA == 'sigma':
                        xlabel = r"$\xi$"
                    else:
                        raise Exception()
                plt.xlabel(xlabel, fontsize=label_size)

                if ylabel is None:
                    if PTA == 'xi':
                        ylabel = r"$\sigma$"
                    elif PTA == 'eta':
                        ylabel = r"$\sigma$"
                    elif PTA == 'sigma':
                        ylabel = r"$\eta$"
                    else:
                        raise Exception()
                plt.ylabel(ylabel, fontsize=label_size)


                plt.tick_params(which='both', labeltop=False, labelright=False, top=True, right=True)
                plt.tick_params(axis='both', which='minor', direction='out', length=minor_tick_length)
                plt.tick_params(axis='both', which='major', direction='out', length=major_tick_length)
                plt.tick_params(axis='both', which='both', labelsize=tick_size)
                plt.tick_params(axis='x', which='both', pad=tick_pad)
                plt.tick_params(axis='y', which='both', pad=tick_pad)


                plt.title(title[n], fontsize=title_size, pad=title_pad)


                if saveto is not None and saveto != '':
                    if self._sf_.k in (0, 3):
                        plt.savefig(saveto, bbox_inches='tight')
                    else:
                        assert saveto.count('.') == 1, f'filename {saveto} is wrong, cannot save to it.'
                        filename, extension = saveto.split('.')
                        plt.savefig(filename + f'_{n}th_component'+'.'+extension, bbox_inches='tight')
                else:
                    plt.show()
                plt.close()
                FIGURES.append(fig)

            return FIGURES, XYZ, VAL

        else:
            return

    def perpendicular_slice_modulus_difference_with(self, MPS, target, normalized=0,
        plot_type='contourf', usetex=True, colormap='coolwarm',
        numOfSamples=100000, figsize=(6, 5),
        num_of_levels=20,
        xlabel=None, ylabel=None,
        title=None, levels=None, # if provide them, put them in list of length 1 (for 0-, 3-form) or 3 (for 1-, 2-form)
        colorbar_label_size=12, colorbar_orientation='vertical', colorbar_ticks=None,
        colorbar_name=None, colorbar_name_pad=0, colorbar_name_size=15,
        title_pad=10,
        label_size=12, title_size=12,
        minor_tick_length=4, major_tick_length=7, tick_pad=8, tick_size=12,
        saveto=None,):
        """If not normalized, let A be a vector (one or two form), A = (u, v, w), and B be another vector, B = (a, b, c), if we do ::

            A.visualize.matplot.perpendicular_slice_sum_square_difference_with(MPS, B)

        where MPS is a _3dCSCG_MeshPerpendicularSlice instance, we actually plot local
            modules_diff = sqrt((u-a)^2 + (v-b)^2 + (w-c)^2)
        over the slice.

        While if A, B are scalars, we then plot local
            modules_diff = abs(A-B).

        The average is defined as, for vectors,
            modules_aver = sqrt(((u+a)/2)^2 + ((v+b)/2)^2 + ((w+c)/2)^2),
        and for scalars,
            modules_aver = abs((A+B)/2)


        :param MPS: A _3dCSCG_MeshPerpendicularSlice instance.
        :param target: compute the difference with ``target``.
        :param normalized:

            If normalized == 0, we plot modules_diff
            If normalized == 1, we divide modules_diff by modules_aver; we plot modules_diff/modules_aver.


        :param plot_type:
        :param figsize:
        :param numOfSamples:
        :param usetex:
        :param colormap:
        :param num_of_levels: Only valid when ``levels`` is ``None``.
        :param xlabel:
        :param ylabel:
        :param title:
        :param levels:
        :param colorbar_label_size:
        :param colorbar_orientation:
        :param colorbar_ticks:
        :param colorbar_name:
        :param colorbar_name_pad:
        :param colorbar_name_size:
        :param title_pad:
        :param label_size:
        :param title_size:
        :param minor_tick_length:
        :param major_tick_length:
        :param tick_pad:
        :param tick_size:
        :param saveto:
        :return:
        """
        assert MPS.__class__.__name__ == "_3dCSCG_MeshPerpendicularSlice", "I need a _3dCSCG_MeshPerpendicularSlice."
        assert MPS._mesh_ == self._sf_.mesh, "Meshes do not match."
        assert MPS._mesh_ == target.mesh, "Meshes do not match."
        assert self._sf_.cochain._local_ is not None, "Form has no cochain!"

        if self._sf_.k in (0, 3):
            assert target.k in (0, 3)
        elif self._sf_.k in (1, 2):
            assert target.k in (1, 2)
        else:
            raise Exception()

        PTA = MPS.perpendicular_to_axis

        loc_len = len(MPS)
        ALL_LEN = COMM.gather(loc_len, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ALL_LEN = sum(ALL_LEN)
            assert ALL_LEN > 0, "MPS is empty!"
            density = int(np.sqrt(numOfSamples / ALL_LEN)) + 1
        else:
            density = None
        density = COMM.bcast(density, root=MASTER_RANK)
        sample = np.linspace(-1, 1, density)

        XYZ, VAL = dict(), dict() # get data for plot.
        for e in MPS:
            eps = MPS[e]
            pta = eps.perpendicular_to_axis
            assert pta == PTA, "For _3dCSCG_MeshPerpendicularSlice, all element perpendicular slice must have same PTA."
            pos = eps.position
            if pta == 'xi':
                xi = np.array([pos,])
                eta = sigma = sample
            elif pta == 'eta':
                eta = np.array([pos,])
                xi = sigma = sample
            elif pta == 'sigma':
                sigma = np.array([pos,])
                xi = eta = sample
            else:
                raise Exception()

            xyz, val = self._sf_.reconstruct(xi, eta, sigma, ravel=False, i=eps._element_.i)

            ___, _a_ = target.reconstruct(xi, eta, sigma, ravel=False, i=eps._element_.i)

            val = val[e]
            xyz = xyz[e]
            _a_ = _a_[e]

            if self._sf_.k in (0, 3):
                assert len(val) == 1
                assert len(_a_) == 1
            else:
                assert len(val) == 3
                assert len(_a_) == 3

            _val_ = list()
            _xyz_ = list()
            _norm = list()

            for i, xyz_i in enumerate(xyz):
                if pta == 'xi':
                    _xyz_.append(xyz_i[0,:,:])
                elif pta == 'eta':
                    _xyz_.append(xyz_i[:,0,:])
                elif pta == 'sigma':
                    _xyz_.append(xyz_i[:,:,0])
                else:
                    raise Exception()

            for i, vi in enumerate(val):

                v_ = _a_[i]

                if pta == 'xi':
                    res = vi[0,:,:]
                    dif = v_[0,:,:]

                elif pta == 'eta':
                    res = vi[:,0,:]
                    dif = v_[:,0,:]

                elif pta == 'sigma':
                    res = vi[:,:,0]
                    dif = v_[:,:,0]

                else:
                    raise Exception()

                _0_ = (res - dif)**2
                _val_.append(_0_)
                _n_ = ((res + dif)/2)**2
                _norm.append(_n_)

            _val_ = np.sqrt(sum(_val_)) # modules of the difference
            if normalized == 0:
                pass
            elif normalized == 1:
                _norm = np.sqrt(sum(_norm)) # modules of the average
                _val_ = _val_/ _norm
            else:
                raise Exception(f"normalized={normalized} is wrong, should be in (0, 1).")

            XYZ[e] = _xyz_
            VAL[e] = _val_

        # gather all information to the mAster core ----------- BELOW ---------------------------------------------
        XYZ = COMM.gather(XYZ, root=MASTER_RANK)
        VAL = COMM.gather(VAL, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ___ = dict()
            for xyz in XYZ:
                ___.update(xyz)
            XYZ = ___

            ___ = dict()
            for val in VAL:
                ___.update(val)
            VAL = ___

            MIN, MAX = None, None

            for e in VAL:
                val = VAL[e]
                min_i = np.min(val)
                max_i = np.max(val)

                if MIN is None:
                    MIN = min_i
                else:
                    if min_i < MIN:
                        MIN = min_i
                    else:
                        pass

                if MAX is None:
                    MAX = max_i
                else:
                    if max_i > MAX:
                        MAX = max_i
                    else:
                        pass

            if MIN == MAX:
                MIN -= 0.5
                MAX += 0.5

            if levels is None:
                levels = np.linspace(MIN, MAX, num_of_levels)


        # Now, we can do the plot ------------- BELOW -----------------------------------------------------------

        if RANK == MASTER_RANK:

            if title is None:

                if self._sf_.k in (0, 3):
                    title = r'$|a-b|$'
                elif self._sf_.k in (1, 2):
                    title = r'$\sqrt{\sum_i(\boldsymbol{u}_i - \boldsymbol{v}_i)^2}$'
                else:
                    raise Exception()

            plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
            plt.rc('text', usetex=usetex)
            if colormap is not None: plt.rcParams['image.cmap'] = colormap


            plotter = getattr(plt, plot_type)

            fig = plt.figure(figsize=figsize)
            # fig, ax = plt.subplots(figsize=figsize)

            for e in VAL: # go through all involved elements.
                if PTA == 'xi':
                    axis_1, axis_2 = XYZ[e][1], XYZ[e][2]
                elif PTA == 'eta':
                    axis_1, axis_2 = XYZ[e][0], XYZ[e][2]
                elif PTA == 'sigma':
                    axis_1, axis_2 = XYZ[e][0], XYZ[e][1]
                else:
                    raise Exception()

                if plot_type == 'contourf':
                    plotter(axis_1, axis_2, VAL[e], levels=levels, extend='both')
                else:
                    plotter(axis_1, axis_2, VAL[e], levels=levels)

            if colorbar_ticks is None: colorbar_ticks = levels
            cb = plt.colorbar(ticks=colorbar_ticks, orientation=colorbar_orientation)
            cb.set_label(colorbar_name, labelpad=colorbar_name_pad, size=colorbar_name_size) # change color bar name, gap and size.
            cb.ax.tick_params(labelsize=colorbar_label_size)
            # cb.ax.set_yticklabels(['< -1', '0', '> 1'])  # We can customized the labels of the colorbar; it must of the same length as ticks.

            if xlabel is None:
                if PTA == 'xi':
                    xlabel = r"$y$"
                elif PTA == 'eta':
                    xlabel = r"$x$"
                elif PTA == 'sigma':
                    xlabel = r"$x$"
                else:
                    raise Exception()
            plt.xlabel(xlabel, fontsize=label_size)

            if ylabel is None:
                if PTA == 'xi':
                    ylabel = r"$z$"
                elif PTA == 'eta':
                    ylabel = r"$z$"
                elif PTA == 'sigma':
                    ylabel = r"$y$"
                else:
                    raise Exception()
            plt.ylabel(ylabel, fontsize=label_size)

            plt.tick_params(which='both', labeltop=False, labelright=False, top=True, right=True)
            plt.tick_params(axis='both', which='minor', direction='out', length=minor_tick_length)
            plt.tick_params(axis='both', which='major', direction='out', length=major_tick_length)
            plt.tick_params(axis='both', which='both', labelsize=tick_size)
            plt.tick_params(axis='x', which='both', pad=tick_pad)
            plt.tick_params(axis='y', which='both', pad=tick_pad)

            plt.title(title, fontsize=title_size, pad=title_pad)

            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')
            else:
                plt.show()
            plt.close()

            return fig, XYZ, VAL

        else:
            return

    def perpendicular_slice_modulus_average_with(self, MPS, target,
        plot_type='contourf', usetex=True, colormap='coolwarm',
        numOfSamples=100000, figsize=(6, 5),
        num_of_levels=20,
        xlabel=None, ylabel=None,
        title=None, levels=None, # if provide them, put them in list of length 1 (for 0-, 3-form) or 3 (for 1-, 2-form)
        colorbar_label_size=12, colorbar_orientation='vertical', colorbar_ticks=None,
        colorbar_name=None, colorbar_name_pad=0, colorbar_name_size=15,
        title_pad=10,
        label_size=12, title_size=12,
        minor_tick_length=4, major_tick_length=7, tick_pad=8, tick_size=12,
        saveto=None,):
        """
        The average is defined as, for vectors,
            modules_aver = sqrt(((u+a)/2)^2 + ((v+b)/2)^2 + ((w+c)/2)^2),
        and for scalars,
            modules_aver = abs((A+B)/2)


        :param MPS: A _3dCSCG_MeshPerpendicularSlice instance.
        :param target: compute the average with ``target``.
        :param plot_type:
        :param figsize:
        :param numOfSamples:
        :param usetex:
        :param colormap:
        :param num_of_levels: Only valid when ``levels`` is ``None``.
        :param xlabel:
        :param ylabel:
        :param title:
        :param levels:
        :param colorbar_label_size:
        :param colorbar_orientation:
        :param colorbar_ticks:
        :param colorbar_name:
        :param colorbar_name_pad:
        :param colorbar_name_size:
        :param title_pad:
        :param label_size:
        :param title_size:
        :param minor_tick_length:
        :param major_tick_length:
        :param tick_pad:
        :param tick_size:
        :param saveto:
        :return:
        """
        assert MPS.__class__.__name__ == "_3dCSCG_MeshPerpendicularSlice", "I need a _3dCSCG_MeshPerpendicularSlice."
        assert MPS._mesh_ == self._sf_.mesh, "Meshes do not match."
        assert MPS._mesh_ == target.mesh, "Meshes do not match."
        assert self._sf_.cochain._local_ is not None, "Form has no cochain!"

        if self._sf_.k in (0, 3):
            assert target.k in (0, 3)
        elif self._sf_.k in (1, 2):
            assert target.k in (1, 2)
        else:
            raise Exception()

        PTA = MPS.perpendicular_to_axis

        loc_len = len(MPS)
        ALL_LEN = COMM.gather(loc_len, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ALL_LEN = sum(ALL_LEN)
            assert ALL_LEN > 0, "MPS is empty!"
            density = int(np.sqrt(numOfSamples / ALL_LEN)) + 1
        else:
            density = None
        density = COMM.bcast(density, root=MASTER_RANK)
        sample = np.linspace(-1, 1, density)

        XYZ, VAL = dict(), dict() # get data for plot.
        for e in MPS:
            eps = MPS[e]
            pta = eps.perpendicular_to_axis
            assert pta == PTA, "For _3dCSCG_MeshPerpendicularSlice, all element perpendicular slice must have same PTA."
            pos = eps.position
            if pta == 'xi':
                xi = np.array([pos,])
                eta = sigma = sample
            elif pta == 'eta':
                eta = np.array([pos,])
                xi = sigma = sample
            elif pta == 'sigma':
                sigma = np.array([pos,])
                xi = eta = sample
            else:
                raise Exception()

            xyz, val = self._sf_.reconstruct(xi, eta, sigma, ravel=False, i=eps._element_.i)

            ___, _a_ = target.reconstruct(xi, eta, sigma, ravel=False, i=eps._element_.i)

            val = val[e]
            xyz = xyz[e]
            _a_ = _a_[e]

            if self._sf_.k in (0, 3):
                assert len(val) == 1
                assert len(_a_) == 1
            else:
                assert len(val) == 3
                assert len(_a_) == 3

            _val_ = list()
            _xyz_ = list()
            _norm = list()

            for i, xyz_i in enumerate(xyz):
                if pta == 'xi':
                    _xyz_.append(xyz_i[0,:,:])
                elif pta == 'eta':
                    _xyz_.append(xyz_i[:,0,:])
                elif pta == 'sigma':
                    _xyz_.append(xyz_i[:,:,0])
                else:
                    raise Exception()

            for i, vi in enumerate(val):

                v_ = _a_[i]

                if pta == 'xi':
                    res = vi[0,:,:]
                    dif = v_[0,:,:]

                elif pta == 'eta':
                    res = vi[:,0,:]
                    dif = v_[:,0,:]

                elif pta == 'sigma':
                    res = vi[:,:,0]
                    dif = v_[:,:,0]

                else:
                    raise Exception()

                _0_ = ((res + dif)/2)**2
                _val_.append(_0_)

            _val_ = np.sqrt(sum(_val_)) # modules of the difference
            XYZ[e] = _xyz_
            VAL[e] = _val_

        # gather all information to the mAster core ----------- BELOW ---------------------------------------------
        XYZ = COMM.gather(XYZ, root=MASTER_RANK)
        VAL = COMM.gather(VAL, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            ___ = dict()
            for xyz in XYZ:
                ___.update(xyz)
            XYZ = ___

            ___ = dict()
            for val in VAL:
                ___.update(val)
            VAL = ___

            MIN, MAX = None, None

            for e in VAL:
                val = VAL[e]
                min_i = np.min(val)
                max_i = np.max(val)

                if MIN is None:
                    MIN = min_i
                else:
                    if min_i < MIN:
                        MIN = min_i
                    else:
                        pass

                if MAX is None:
                    MAX = max_i
                else:
                    if max_i > MAX:
                        MAX = max_i
                    else:
                        pass

            if MIN == MAX:
                MIN -= 0.5
                MAX += 0.5

            if levels is None:
                levels = np.linspace(MIN, MAX, num_of_levels)


        # Now, we can do the plot ------------- BELOW -----------------------------------------------------------

        if RANK == MASTER_RANK:

            if title is None:

                if self._sf_.k in (0, 3):
                    title = r'$|(a+b)/2|$'
                elif self._sf_.k in (1, 2):
                    title = r'$\left\|\dfrac{\boldsymbol{u}+\boldsymbol{v}}{2}\right\|$'
                else:
                    raise Exception()

            plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
            plt.rc('text', usetex=usetex)
            if colormap is not None: plt.rcParams['image.cmap'] = colormap


            plotter = getattr(plt, plot_type)

            fig = plt.figure(figsize=figsize)
            # fig, ax = plt.subplots(figsize=figsize)

            for e in VAL: # go through all involved elements.
                if PTA == 'xi':
                    axis_1, axis_2 = XYZ[e][1], XYZ[e][2]
                elif PTA == 'eta':
                    axis_1, axis_2 = XYZ[e][0], XYZ[e][2]
                elif PTA == 'sigma':
                    axis_1, axis_2 = XYZ[e][0], XYZ[e][1]
                else:
                    raise Exception()

                if plot_type == 'contourf':
                    plotter(axis_1, axis_2, VAL[e], levels=levels, extend='both')
                else:
                    plotter(axis_1, axis_2, VAL[e], levels=levels)

            if colorbar_ticks is None: colorbar_ticks = levels
            cb = plt.colorbar(ticks=colorbar_ticks, orientation=colorbar_orientation)
            cb.set_label(colorbar_name, labelpad=colorbar_name_pad, size=colorbar_name_size) # change color bar name, gap and size.
            cb.ax.tick_params(labelsize=colorbar_label_size)
            # cb.ax.set_yticklabels(['< -1', '0', '> 1'])  # We can customized the labels of the colorbar; it must of the same length as ticks.

            if xlabel is None:
                if PTA == 'xi':
                    xlabel = r"$y$"
                elif PTA == 'eta':
                    xlabel = r"$x$"
                elif PTA == 'sigma':
                    xlabel = r"$x$"
                else:
                    raise Exception()
            plt.xlabel(xlabel, fontsize=label_size)

            if ylabel is None:
                if PTA == 'xi':
                    ylabel = r"$z$"
                elif PTA == 'eta':
                    ylabel = r"$z$"
                elif PTA == 'sigma':
                    ylabel = r"$y$"
                else:
                    raise Exception()
            plt.ylabel(ylabel, fontsize=label_size)

            plt.tick_params(which='both', labeltop=False, labelright=False, top=True, right=True)
            plt.tick_params(axis='both', which='minor', direction='out', length=minor_tick_length)
            plt.tick_params(axis='both', which='major', direction='out', length=major_tick_length)
            plt.tick_params(axis='both', which='both', labelsize=tick_size)
            plt.tick_params(axis='x', which='both', pad=tick_pad)
            plt.tick_params(axis='y', which='both', pad=tick_pad)

            plt.title(title, fontsize=title_size, pad=title_pad)

            if saveto is not None and saveto != '':
                plt.savefig(saveto, bbox_inches='tight')
            else:
                plt.show()
            plt.close()

            return fig, XYZ, VAL

        else:
            return



if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\standard\visualize\matplot.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([2,3,4])
    space = SpaceInvoker('polynomials')([('Lobatto',4), ('Lobatto',3), ('Lobatto',2)])
    FC = FormCaller(mesh, space)

    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    f0 = FC('0-f', is_hybrid=False, numbering_parameters={'scheme_name': 'Naive',})
    f0.TW.func.do.set_func_body_as(es, 'pressure')
    f0.TW.current_time = 0
    f0.TW.do.push_all_to_instant()
    f0.do.discretize()

    f1 = FC('1-f', is_hybrid=False, numbering_parameters={'scheme_name': 'Naive',})
    f1.TW.func.do.set_func_body_as(es, 'velocity')
    f1.TW.current_time = 0
    f1.TW.do.push_all_to_instant()
    f1.do.discretize()

    f2 = FC('2-f', is_hybrid=False, numbering_parameters={'scheme_name': 'Naive',})
    f2.TW.func.do.set_func_body_as(es, 'velocity')
    f2.TW.current_time = 0
    f2.TW.do.push_all_to_instant()
    f2.do.discretize()

    f3 = FC('3-f', is_hybrid=False, numbering_parameters={'scheme_name': 'Naive',})
    f3.TW.func.do.set_func_body_as(es, 'pressure')
    f3.TW.current_time = 0
    f3.TW.do.push_all_to_instant()
    f3.do.discretize()

    region = mesh.domain.regions['R:R']
    RS = region.sub_geometry.make_a_perpendicular_slice_object_on(r=4.5 / 9)
    MS = mesh.sub_geometry.make_a_perpendicular_slice_object_on(RS)

    f0.visualize.matplot.perpendicular_slice(MS)

    # f0.visualize.matplot.perpendicular_slice_sum_square_difference_with(MS, f3, saveto='')
    # f1.visualize.matplot.perpendicular_slice_modules_difference_with(MS, f2, saveto='', normalized=0)
    # f1.visualize.matplot.perpendicular_slice_modules_average_with(MS, f2, saveto='')
    # f2.visualize.matplot.perpendicular_slice_sum_square_difference_with(MS, f1, saveto='')