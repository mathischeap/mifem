# -*- coding: utf-8 -*-

from objects.CSCG._3d.forms.standard.base.visualize.matplot import _3dCSCG_standard_form_Matplot
from root.config.main import *
import matplotlib.pyplot as plt
import matplotlib


class _3dCSCG_S3F_VISUALIZE_Matplot(_3dCSCG_standard_form_Matplot):
    """"""
    def __init__(self, sf):
        """"""
        super(_3dCSCG_S3F_VISUALIZE_Matplot, self).__init__(sf)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        return self.perpendicular_surface(*args, **kwargs)

    def perpendicular_surface(self,
        x=None, y=None, z=None, # only one of them can be not-None.
        plot_type='contourf', usetex=True, colormap='coolwarm',
        numOfSamples=100000, figsize=(6, 5),
        num_of_levels=20,
        xlabel=None, ylabel=None,
        title=None, levels=None, # if provide them, put them in list of length 1 (for 0-, 3-form) or 3 (for 1-, 2-form)
        colorbar_font_size=12, title_pad=10,
        label_size=12, title_size=12,
        minor_tick_length=4, major_tick_length=7, tick_pad=8, tick_size=12,
        saveto=None,):
        """Plot the 3-sf on a surface of x=constant or y=constant or z=constant.

        Note that the plot will only be reasonable if the region mappings are regular. For example,
        see `region.interpolation.___inverse_mapping_r_x_s0t0___`.

        Parameters
        ----------
        x
        y
        z
        plot_type
        usetex
        colormap
        numOfSamples
        figsize
        num_of_levels
        xlabel
        ylabel
        title
        levels
        colorbar_font_size
        title_pad
        label_size
        title_size
        minor_tick_length
        major_tick_length
        tick_pad
        tick_size
        saveto

        Returns
        -------

        """
        assert self._sf_.cochain._local_ is not None, "Form has no cochain!"

        s3f = self._sf_
        mesh = s3f.mesh
        MPS = mesh.sub_geometry.make_a_perpendicular_slice_object_on(x=x, y=y, z=z)
        assert MPS.__class__.__name__ == "_3dCSCG_MeshPerpendicularSlice"

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

            xyz, val = self._sf_.reconstruct(xi, eta, sigma, ravel=False, element_range=eps._element_.i)
            val = val[e]
            xyz = xyz[e]

            assert len(val) == 1

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

            MIN, MAX = None, None

            for e in VAL:
                for val in VAL[e]:
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
                levels= np.linspace(MIN, MAX, num_of_levels)


        # Now, we can do the plot ------------- BELOW -----------------------------------------------------------

        if RANK == MASTER_RANK:

            if title is None:
                title = f'{self._sf_.k}-form: {self._sf_.standard_properties.name}'

            if saveto is not None: matplotlib.use('Agg')

            plt.rc('text', usetex=usetex)
            if colormap is not None: plt.rcParams['image.cmap'] = colormap

            plotter = getattr(plt, plot_type)

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

                plotter(axis_1, axis_2, VAL[e][0], levels=levels)

            cb = plt.colorbar()
            # cb.set_label("Relative Photon Intensity", labelpad=-1, size=14) # change color bar name, gap and size.
            cb.ax.tick_params(labelsize=colorbar_font_size)

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
                pass

            plt.show()
            plt.close()

            return fig, XYZ, VAL

        else:
            return