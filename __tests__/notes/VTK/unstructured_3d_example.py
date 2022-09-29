
import os
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkVoxel, VtkHexahedron # the 11th and 12th type of unstructured VTK cell.
import numpy as np

FILE_PATH = "./unstructured"


def clean():
    try:
        os.remove(FILE_PATH + ".vtu")
    except:
        pass


def run():

    # Define vertices
    x = np.zeros(12)
    y = np.zeros(12)
    z = np.zeros(12)

    #
    x[0], y[0], z[0] = 0.0, 0.0, 0.0
    x[1], y[1], z[1] = 1.0, 0.0, 0.0
    x[2], y[2], z[2] = 1.0, 1.0, 0.0
    x[3], y[3], z[3] = 0.0, 1.0, 0.0
    x[4], y[4], z[4] = 0.0, 0.0, 1.0
    x[5], y[5], z[5] = 1.0, 0.0, 1.0
    x[6], y[6], z[6] = 1.0, 1.0, 1.0
    x[7], y[7], z[7] = 0.0, 1.0, 1.0
    x[8], y[8], z[8] = 1.1, 2.0, 0.2
    x[9], y[9], z[9] = 0.1, 2.0, 0.2
    x[10], y[10], z[10] = 1.1, 2.0, 1.2
    x[11], y[11], z[11] = 0.1, 2.0, 1.2

    # Define connectivity or vertices that belongs to each element
    conn = np.zeros(16)
    #
    conn[:8] = [0,1,2,3,4,5,6,7]  # first hexahedron
    conn[8:] = [3,2,8,9,7,6,10,11]  # second hexahedron
    #
    # Define offset of last vertex of each element
    offset = np.zeros(2)
    offset[0] = 8
    offset[1] = 16

    # Define cell types

    ctype = np.zeros(2)
    ctype[0], ctype[1] = VtkHexahedron.tid, VtkHexahedron.tid
    # ctype[2] = VtkQuad.tid
    #
    # cd = np.random.rand(3)
    # cellData = {"pressure": cd}
    #
    pd = np.random.rand(12)
    pointData = {"temp": pd}
    #
    print(x, y, z)
    print(conn)
    print(offset)
    print(ctype)
    unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity=conn, offsets=offset, cell_types=ctype,
                          # cellData=cellData,
                          pointData=pointData)


if __name__ == "__main__":
    # python __tests__/notes/VTK/unstructured_3d_example.py
    clean()
    run()