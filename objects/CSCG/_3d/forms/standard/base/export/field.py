# -*- coding: utf-8 -*-
"""We want to export the field to some data files.
"""

from root.config.main import *
from screws.freeze.main import FrozenOnly
from screws.miscellaneous.timer import check_filename, check_no_splcharacter
from scipy.io import savemat


class _3dCSC_SF_Export_Field(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        assert '3dCSCG_standard_form' in sf.standard_properties.tags
        self._sf_ = sf
        self._freeze_self_()


    def to_file(self, filename, numOfSamples=1e6, regions=None):
        """"""
        filename, extension = check_filename(filename)
        if extension is None: extension = 'txt'

        supported_formats = ('txt', 'mat')
        assert extension in supported_formats, \
            f"format={extension} is not among the supported formats {supported_formats}."

        if isinstance(numOfSamples, (int, float)):
            assert numOfSamples > 0, f"numOfSamples={numOfSamples} is wrong."
            numOfSamples = [numOfSamples, numOfSamples, numOfSamples]
        else:
            assert isinstance(numOfSamples, (tuple, list)) and len(numOfSamples) == 3, \
                f"numOfSamples={numOfSamples} wrong."
            for nos in numOfSamples:
                assert isinstance(nos, (int, float)) and nos > 0, f"numOfSamples={numOfSamples} wrong."

        mesh = self._sf_.mesh

        if regions is None:
            regions = mesh.domain.regions.names
        elif isinstance(regions, str):
            regions = [regions,]
        else:
            pass
        assert isinstance(regions, (list, tuple)), f"regions={regions} is wrong."
        assert len(set(regions)) == len(regions), f"regions={regions} has repeated regions."
        for i, r in enumerate(regions):
            assert r in mesh.domain.regions, f"regions[{i}]={r} is wrong."

        rst = list()
        for i in range(3):
            density = int((numOfSamples[i] / mesh.elements.GLOBAL_num) ** (1/3)) + 1
            interval = 2 / density
            rst.append(np.linspace(-1 + interval/2, 1-interval/2, density))

        xyz, v = self._sf_.reconstruct(*rst, regions=regions)

        # Now, we gather xyz & v from all cores into Master Core, store in XYZ & V --- BELOW ---
        if RANK == MASTER_RANK:
            X = [None for _ in range(mesh.elements.GLOBAL_num)]
            Y = [None for _ in range(mesh.elements.GLOBAL_num)]
            Z = [None for _ in range(mesh.elements.GLOBAL_num)]
            Vx = [None for _ in range(mesh.elements.GLOBAL_num)]
            if self._sf_.k in (1, 2):
                Vy = [None for _ in range(mesh.elements.GLOBAL_num)]
                Vz = [None for _ in range(mesh.elements.GLOBAL_num)]
            for j in mesh.elements.indices:
                X[j] = xyz[j][0]
                Y[j] = xyz[j][1]
                Z[j] = xyz[j][2]
                Vx[j] = v[j][0]
                if self._sf_.k in (1, 2):
                    # noinspection PyUnboundLocalVariable
                    Vy[j] = v[j][1]
                    # noinspection PyUnboundLocalVariable
                    Vz[j] = v[j][2]
            for i in SLAVE_RANKS:
                xyz, v = COMM.recv(source=i, tag=0)
                for j in xyz:
                    X[j] = xyz[j][0]
                    Y[j] = xyz[j][1]
                    Z[j] = xyz[j][2]
                    Vx[j] = v[j][0]
                    if self._sf_.k in (1, 2):
                        Vy[j] = v[j][1]
                        Vz[j] = v[j][2]
            del xyz, v
        else:
            COMM.send([xyz, v], dest=MASTER_RANK, tag=0)
            del xyz, v

        # Now, we reshape the XYZ and V for export in the master core. -------- BELOW ----------
        if RANK == MASTER_RANK:
            if self._sf_.k in (1, 2):
                # noinspection PyUnboundLocalVariable
                X, Y, Z, Vx, Vy, Vz = mesh.do.regionwsie_stack(X, Y, Z, Vx, Vy, Vz)
            else:
                # noinspection PyUnboundLocalVariable
                X, Y, Z, V = mesh.do.regionwsie_stack(X, Y, Z, Vx)

            for rn in regions:
                assert rn in X and rn in Y and rn in Z, "Data not full!"

                x, y, z = X[rn], Y[rn], Z[rn]
                if self._sf_.k in (1, 2):
                    vx, vy, vz = Vx[rn], Vy[rn], Vz[rn]
                else:
                    # noinspection PyUnboundLocalVariable
                    vx = V[rn]

                # we take care of the file names ------------------ BELOW -----------------------
                RN = rn[2:] # if regions name is R:center, we select
                assert check_no_splcharacter(RN), f"region name={RN} wrong."

                FILE_NAME = filename + '__InRegion_' + RN
                if self._sf_.k in (1, 2):
                    FILE_NAME += '__x_y_z_vx_vy_vz'
                else:
                    FILE_NAME += '__x_y_z_v'
                FILE_NAME = FILE_NAME + '.' + extension


                # It's time to do the save or writing ------------------- BELOW -----------------

                if extension == 'txt':
                    # for .txt, we have to flat the data =====================
                    x = x.ravel(order='F')[:,np.newaxis]
                    y = y.ravel(order='F')[:,np.newaxis]
                    z = z.ravel(order='F')[:,np.newaxis]
                    if self._sf_.k in (1, 2):
                        vx = vx.ravel(order='F')[:,np.newaxis]
                        # noinspection PyUnboundLocalVariable
                        vy = vy.ravel(order='F')[:,np.newaxis]
                        # noinspection PyUnboundLocalVariable
                        vz = vz.ravel(order='F')[:,np.newaxis]
                    else:
                        vx = vx.ravel(order='F')[:,np.newaxis]
                    if self._sf_.k in (1, 2):
                        # noinspection PyUnboundLocalVariable
                        TO_BE_WRITTEN = np.hstack((x, y, z, vx, vy, vz))
                    else:
                        TO_BE_WRITTEN = np.hstack((x, y, z, vx))
                    # noinspection PyTypeChecker
                    np.savetxt(FILE_NAME, TO_BE_WRITTEN)

                elif extension == 'mat':
                    # for .mat, we save 3-d arrays. ==========================
                    m_dic = dict()
                    m_dic['x'] = x
                    m_dic['y'] = y
                    m_dic['z'] = z
                    if self._sf_.k in (1, 2):
                        m_dic['vx'] = vx
                        m_dic['vy'] = vy
                        m_dic['vz'] = vz
                    else:
                        m_dic['v'] = vx

                    savemat(FILE_NAME, m_dic)

                else:
                    raise Exception(f"Format=.{extension} is not supported.")