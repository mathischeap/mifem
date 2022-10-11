# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/02 5:25 PM
"""
import sys
import numpy as np

from evtk import hl
# from pyevtk.vtk import VtkVoxel, VtkHexahedron # the 11th and 12th type of unstructured VTK cell.
from scipy.sparse import lil_matrix
if './' not in sys.path: sys.path.append('./')

from root.config.main import cOmm, rAnk, mAster_rank

from objects.CSCG._3d.forms.standard._0s.main import _3dCSCG_0Form
from objects.CSCG._3d.spaces.polynomials import _3dCSCG_PolynomialSpace


from objects.CSCG._2d.forms.standard._0_form.outer.main import _2dCSCG_0Form_Outer
from objects.CSCG._2d.spaces.polynomials import _2dCSCG_PolynomialSpace


def unstructuredGridToVTK(grid, objs, filename):
    """This works for all types of domains. But it is slower. For domain of a single (or a few)
    region(s), you may want to use `gridToVTK` instead.

    This function is ideal for domains of multiple regions. For domains of only one single (ro
    a few) regions(s), it is recommended to use `gridToVTK`.

    Parameters
    ----------
    grid : {list, dict}
        If it is a list, we will make it into a dict whose keys are the region names  of the mesh,
        and this list is the value for all keys.

        For each value, if the mesh is 3 (2) d, it is a list of three (two) 1d array.

        The mesh will be the mesh of the `objs`. So we need all objs to be of the same mesh.

        IMPORTANT: please notice that this grid must be conforming across regions!

    objs :
        The CSCG objects to be visualized.
    filename :

    Returns
    -------

    """
    if not isinstance(objs, (list, tuple)):
        objs = [objs,]
    else:
        pass
    for i, obj in enumerate(objs):
        assert hasattr(obj, 'mesh'),  f"{i}th obj does not have a mesh."
        if i == 0:
            mesh = obj.mesh
        else:
            # noinspection PyUnboundLocalVariable
            assert mesh == obj.mesh, f"mesh of {i}th obj does not match that of 0th obj."

    DiscreteFields = list()
    for i, obj in enumerate(objs):
        df = obj.reconstruct.discrete_field(grid)
        if i == 0:
            xyz = df.coordinates
        else:
            _xyz = df.coordinates
            # noinspection PyUnboundLocalVariable
            for rn in xyz:
                for j, axis in enumerate(xyz[rn]):
                    np.testing.assert_array_almost_equal(axis, _xyz[rn][j], decimal=8)
            df._coordinates_ = None # clean COO as we only need it in df[0]
        DiscreteFields.append(df)
    del _xyz

    ndim = mesh.ndim
    if ndim == 3:
        return _3dCSCG_unstructuredGridToVTK(mesh, DiscreteFields, filename, objs)
    elif ndim == 2:
        return _2dCSCG_unstructuredGridToVTK(mesh, DiscreteFields, filename, objs)
    else:
        raise Exception





def _2dCSCG_unstructuredGridToVTK(mesh, dfs, filename, objs):
    """The function for objects in 2d CSCG meshes.

    Parameters
    ----------
    mesh
    dfs
    filename
    objs

    Returns
    -------

    """
    #------ get element_layout for the new reference mesh -----------------------------------------1
    df0 = dfs[0]
    GRID = cOmm.gather(df0.grid, root=mAster_rank)
    if rAnk == mAster_rank:
        _grid = dict()
        for G in GRID:
            _grid.update(G)
        GRID = _grid
    GRID = cOmm.bcast(GRID, root=mAster_rank)
    element_layout = dict()
    for rn in GRID:
        element_layout[rn] = [
            np.diff(_) for _ in GRID[rn]
        ]
    new_mesh = mesh.__class__(mesh.domain, element_layout=element_layout) #-- produce new mesh ----1

    new_space = _2dCSCG_PolynomialSpace([1,1], None)   #-------------------- produce new space ----1

    r0f = _2dCSCG_0Form_Outer(new_mesh, new_space, is_hybrid=False) # produce the reference 0-form 1

    GM = r0f.numbering.gathering #----- the gathering numbering, a key property we want to use ----1
    Gnd = GM.GLOBAL_num_dofs # -------- this many 0-form dofs = this many nodes in the VTK --------1

    total_num_cells = new_mesh.elements.GLOBAL_num

    #----- collect needed COO ---------------------------------------------------------------------1
    COO = df0.coordinates
    in_regions = new_mesh.elements.in_regions

    ALL_COO = cOmm.allgather(COO)
    need_COO = dict()
    for ac in ALL_COO:
        for bn in ac:
            if bn in in_regions:
                need_COO[bn] = ac[bn]
    COO = need_COO
    del ALL_COO
    df0._coordinates_ = None # clean

    #------ variable names ------------------------------------------------------------------------1
    var_names = list()
    for i, obj in enumerate(objs):
        if hasattr(obj, 'name'):
            name = obj.name
        elif hasattr(obj, 'standard_properties'):
            name = str(obj.standard_properties.name)
        else:
            name = f'variable_{i}'

        j = 0
        old_name = name
        while name in var_names:
            j += 1
            name = old_name + '_' + str(j)

        var_names.append(name)

    #------ collect values ------------------------------------------------------------------------1

    var_dim = dict()
    ALL_needed_VAL = dict()
    for i, df in enumerate(dfs):
        VAL = df.values

        VAL = cOmm.allgather(VAL)
        need_VAL = dict()
        for v in VAL:
            for bn in v:
                if bn in in_regions:
                    need_VAL[bn] = v[bn]
        del VAL
        df._values_ = None
        name = var_names[i]
        ALL_needed_VAL[name] = need_VAL
        var_dim[name] = df.vdim

    #----------------------------- lets make the VTK ----------------------------------------------1
    X, Y = lil_matrix((1, Gnd)), lil_matrix((1, Gnd))

    VAR = dict()
    for vn in var_names:
        VAR[vn] = list()
        for _ in range(var_dim[vn]):
            VAR[vn].append(lil_matrix((1, Gnd)))

    CONN = lil_matrix((total_num_cells, 4), dtype=int)
    cTypes = lil_matrix((1, total_num_cells), dtype=int)

    swift_matrix = np.array([
        [1, 0, 0, 0], # 0 <- 0
        [0, 1, 0, 0], # 1 <- 1
        [0, 0, 0, 1], # 2 <- 3
        [0, 0, 1, 0], # 3 <- 2
    ])

    for e in new_mesh.elements:
        element = new_mesh.elements[e]
        rn, indices = new_mesh.do.find.region_name_and_local_indices_of_element(e)
        x, y, = COO[rn]
        gv = GM[e].full_vector
        i, j = indices

        I_seq = [i, i+1, i, i+1]
        J_seq = [j, j, j+1, j+1]

        X[0, gv] = x[I_seq, J_seq]
        Y[0, gv] = y[I_seq, J_seq]

        for vn in var_names:
            v = ALL_needed_VAL[vn][rn]
            for _ in range(var_dim[vn]):
                VAR[vn][_][0, gv] = v[_][I_seq, J_seq]

        if element.IS.orthogonal:
            CONN[e, :] = gv
            cTypes[0, e] = 8 # VTK PIXEL
        else:
            CONN[e, :] = swift_matrix @ gv
            cTypes[0, e] = 9 # VTK _QUAD

    del ALL_needed_VAL, COO

    X = cOmm.gather(X, root=mAster_rank)
    Y = cOmm.gather(Y, root=mAster_rank)

    VAR = cOmm.gather(VAR, root=mAster_rank)

    CONN = cOmm.gather(CONN, root=mAster_rank)
    cTypes = cOmm.gather(cTypes, root=mAster_rank)

    if rAnk == mAster_rank:
        _x, _y = np.zeros(Gnd), np.zeros(Gnd)
        for x_ in X:
            x_ = x_.tocsr()
            data = x_.data
            indices = x_.indices
            _x[indices] = data
        for y_ in Y:
            y_ = y_.tocsr()
            data = y_.data
            indices = y_.indices
            _y[indices] = data

        V = dict()
        for vn in var_names:
            V[vn] = list()
            for _ in range(var_dim[vn]):
                V[vn].append(np.zeros(Gnd))

        for var_core in VAR:
            for vn in var_core:
                for _ in range(var_dim[vn]):
                    v_ = var_core[vn][_].tocsr()
                    data = v_.data
                    indices = v_.indices
                    V[vn][_][indices] = data

        # noinspection PyUnresolvedReferences
        CONN  = np.sum(CONN, axis=0).toarray().ravel('C')
        offsets = np.linspace(4, CONN.size, total_num_cells)

        # noinspection PyUnresolvedReferences
        cTypes  = np.sum(cTypes, axis=0).toarray().ravel()

        for vn in V:
            if  len(V[vn]) == 1:
                # noinspection PyUnresolvedReferences
                V[vn] = V[vn][0]
            else:
                V[vn].append(np.zeros_like(V[vn][0]))
                # noinspection PyUnresolvedReferences
                V[vn] = tuple(V[vn]) # do not use list

        full_save_path = hl.unstructuredGridToVTK(filename, _x, _y, np.zeros_like(_x),
                                 connectivity=CONN, offsets=offsets, cell_types=cTypes,
                                 pointData=V
                                 )

    else:
        full_save_path = None

    return cOmm.bcast(full_save_path, root=mAster_rank)





def _3dCSCG_unstructuredGridToVTK(mesh, dfs, filename, objs):
    """The function for objects in 3d CSCG meshes.

    Parameters
    ----------
    mesh
    dfs
    filename
    objs

    Returns
    -------

    """
    #------ get element_layout for the new reference mesh -----------------------------------------1
    df0 = dfs[0]
    GRID = cOmm.gather(df0.grid, root=mAster_rank)
    if rAnk == mAster_rank:
        _grid = dict()
        for G in GRID:
            _grid.update(G)
        GRID = _grid
    GRID = cOmm.bcast(GRID, root=mAster_rank)
    element_layout = dict()
    for rn in GRID:
        element_layout[rn] = [
            np.diff(_) for _ in GRID[rn]
        ]
    new_mesh = mesh.__class__(mesh.domain, element_layout=element_layout) #-- produce new mesh ----1

    new_space = _3dCSCG_PolynomialSpace([1,1,1], None) #-------------------- produce new space ----1

    r0f = _3dCSCG_0Form(new_mesh, new_space, is_hybrid=False) #-- produce the reference 0-form ----1

    GM = r0f.numbering.gathering #----- the gathering numbering, a key property we want to use ----1
    Gnd = GM.GLOBAL_num_dofs # -------- this many 0-form dofs = this many nodes in the VTK --------1

    total_num_cells = new_mesh.elements.GLOBAL_num

    #----- collect needed COO ---------------------------------------------------------------------1
    COO = df0.coordinates
    in_regions = new_mesh.elements.in_regions

    ALL_COO = cOmm.allgather(COO)
    need_COO = dict()
    for ac in ALL_COO:
        for bn in ac:
            if bn in in_regions:
                need_COO[bn] = ac[bn]
    COO = need_COO
    del ALL_COO
    df0._coordinates_ = None # clean

    #------ variable names ------------------------------------------------------------------------1
    var_names = list()
    for i, obj in enumerate(objs):
        if hasattr(obj, 'name'):
            name = obj.name
        elif hasattr(obj, 'standard_properties'):
            name = str(obj.standard_properties.name)
        else:
            name = f'variable_{i}'

        j = 0
        old_name = name
        while name in var_names:
            j += 1
            name = old_name + '_' + str(j)

        var_names.append(name)

    #------ collect values ------------------------------------------------------------------------1

    var_dim = dict()
    ALL_needed_VAL = dict()
    for i, df in enumerate(dfs):
        VAL = df.values

        VAL = cOmm.allgather(VAL)
        need_VAL = dict()
        for v in VAL:
            for bn in v:
                if bn in in_regions:
                    need_VAL[bn] = v[bn]
        del VAL
        df._values_ = None
        name = var_names[i]
        ALL_needed_VAL[name] = need_VAL
        var_dim[name] = df.vdim

    #----------------------------- lets make the VTK ----------------------------------------------1
    X, Y, Z = lil_matrix((1, Gnd)), lil_matrix((1, Gnd)), lil_matrix((1, Gnd))

    VAR = dict()
    for vn in var_names:
        VAR[vn] = list()
        for _ in range(var_dim[vn]):
            VAR[vn].append(lil_matrix((1, Gnd)))

    CONN = lil_matrix((total_num_cells, 8), dtype=int)
    cTypes = lil_matrix((1, total_num_cells), dtype=int)

    swift_matrix = np.array([
        [1, 0, 0, 0, 0, 0, 0, 0], # 0 <- 0
        [0, 1, 0, 0, 0, 0, 0, 0], # 1 <- 1
        [0, 0, 0, 1, 0, 0, 0, 0], # 2 <- 3
        [0, 0, 1, 0, 0, 0, 0, 0], # 3 <- 2
        [0, 0, 0, 0, 1, 0, 0, 0], # 4 <- 4
        [0, 0, 0, 0, 0, 1, 0, 0], # 5 <- 5
        [0, 0, 0, 0, 0, 0, 0, 1], # 6 <- 7
        [0, 0, 0, 0, 0, 0, 1, 0], # 7 <- 6
    ])

    for e in new_mesh.elements:
        element = new_mesh.elements[e]
        rn, indices = new_mesh.do.find.region_name_and_local_indices_of_element(e)
        x, y, z = COO[rn]
        gv = GM[e].full_vector
        i, j, k = indices

        I_seq = [i, i+1, i, i+1, i, i+1, i, i+1]
        J_seq = [j, j, j+1, j+1, j, j, j+1, j+1]
        K_seq = [k, k, k, k, k+1, k+1, k+1, k+1]

        X[0, gv] = x[I_seq, J_seq, K_seq]
        Y[0, gv] = y[I_seq, J_seq, K_seq]
        Z[0, gv] = z[I_seq, J_seq, K_seq]

        for vn in var_names:
            v = ALL_needed_VAL[vn][rn]
            for _ in range(var_dim[vn]):
                VAR[vn][_][0, gv] = v[_][I_seq, J_seq, K_seq]

        if element.IS.orthogonal:
            CONN[e, :] = gv
            cTypes[0, e] = 11 # VTK Voxel
        else:
            CONN[e, :] = swift_matrix @ gv
            cTypes[0, e] = 12 # VTK Hexahedron

    del ALL_needed_VAL, COO

    X = cOmm.gather(X, root=mAster_rank)
    Y = cOmm.gather(Y, root=mAster_rank)
    Z = cOmm.gather(Z, root=mAster_rank)

    VAR = cOmm.gather(VAR, root=mAster_rank)

    CONN = cOmm.gather(CONN, root=mAster_rank)
    cTypes = cOmm.gather(cTypes, root=mAster_rank)

    if rAnk == mAster_rank:
        _x, _y, _z = np.zeros(Gnd), np.zeros(Gnd), np.zeros(Gnd)
        for x_ in X:
            x_ = x_.tocsr()
            data = x_.data
            indices = x_.indices
            _x[indices] = data
        for y_ in Y:
            y_ = y_.tocsr()
            data = y_.data
            indices = y_.indices
            _y[indices] = data
        for z_ in Z:
            z_ = z_.tocsr()
            data = z_.data
            indices = z_.indices
            _z[indices] = data

        V = dict()
        for vn in var_names:
            V[vn] = list()
            for _ in range(var_dim[vn]):
                # noinspection PyUnresolvedReferences
                V[vn].append(np.zeros(Gnd))

        for var_core in VAR:
            for vn in var_core:
                for _ in range(var_dim[vn]):
                    v_ = var_core[vn][_].tocsr()
                    data = v_.data
                    indices = v_.indices
                    V[vn][_][indices] = data

        # noinspection PyUnresolvedReferences
        CONN  = np.sum(CONN, axis=0).toarray().ravel('C')
        offsets = np.linspace(8, CONN.size, total_num_cells)

        # noinspection PyUnresolvedReferences
        cTypes  = np.sum(cTypes, axis=0).toarray().ravel()

        for vn in V:
            if  len(V[vn]) == 1:
                # noinspection PyUnresolvedReferences
                V[vn] = V[vn][0]
            else:
                # noinspection PyUnresolvedReferences
                V[vn] = tuple(V[vn]) # do not use list
        full_save_path = hl.unstructuredGridToVTK(filename, _x, _y, _z,
                                 connectivity=CONN, offsets=offsets, cell_types=cTypes,
                                 pointData=V
                                 )

    else:
        full_save_path = None

    return cOmm.bcast(full_save_path, root=mAster_rank)



if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/tools/unstructuredGridToVTK.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('bridge_arch_cracked')([8,8,8])
    # mesh = MeshGenerator('crazy', c=0.)([10,10,10])
    # for rn in mesh.domain.regions:
    #     region = mesh.domain.regions[rn]
    #     print(region.IS.orthogonal)

    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    def u(t,x,y,z): return np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y) + t

    def p(t,x,y,z): return np.sin(2*np.pi*x)*np.sin(np.pi*y)*np.sin(2*np.pi*z) + t

    scalar = FC('scalar', p)
    vector = FC('vector', (u,v,w))

    f0 = FC('0-f', is_hybrid=False, name='pressure')
    f0.TW.func.do.set_func_body_as(scalar)
    f0.TW.current_time = 0
    f0.TW.___DO_push_all_to_instant___()
    f0.discretize()

    f1 = FC('1-f', is_hybrid=False, name='vorticity')
    f1.TW.func.do.set_func_body_as(vector)
    f1.TW.current_time = 0
    f1.TW.___DO_push_all_to_instant___()
    f1.discretize()

    f2 = FC('2-f', is_hybrid=False, name='velocity')
    f2.TW.func.do.set_func_body_as(vector)
    f2.TW.current_time = 0
    f2.TW.___DO_push_all_to_instant___()
    f2.discretize()
    f3 = FC('3-f', is_hybrid=False, name='total pressure')
    f3.TW.func.do.set_func_body_as(scalar)
    f3.TW.current_time = 0
    f3.TW.___DO_push_all_to_instant___()
    f3.discretize()

    grid = [np.linspace(-1,1,15), np.linspace(-1,1,15), np.linspace(-1,1,15)]

    unstructuredGridToVTK(grid, [f0, f2], 'unstructuredGridToVTK_test')

    # ----------- 2d test below ------------------------------------------------
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller

    # mesh = MeshGenerator('cic')([3,3])
    # space = SpaceInvoker('polynomials')([3,3])
    # FC = FormCaller(mesh, space)
    #
    # def u(t,x,y): return np.sin(np.pi*x)*np.cos(2*np.pi*y)+t
    # def v(t,x,y): return np.cos(np.pi*x)*np.sin(np.pi*y)+t
    #
    # def p(t,x,y): return np.sin(2*np.pi*x)*np.sin(np.pi*y)+t
    #
    # scalar = FC('scalar', p)
    # velocity = FC('vector', (u,v))
    #
    # f0 = FC('0-f-o', is_hybrid=False, name='pressure')
    # f0.TW.func.do.set_func_body_as(scalar)
    # f0.TW.current_time = 0
    # f0.TW.___DO_push_all_to_instant___()
    # f0.discretize()
    # f1i = FC('1-f-i', is_hybrid=False, name='Velocity')
    # f1i.TW.func.do.set_func_body_as(velocity)
    # f1i.TW.current_time = 0
    # f1i.TW.___DO_push_all_to_instant___()
    # f1i.discretize()
    # f1o = FC('1-f-o', is_hybrid=False, name='vorticity')
    # f1o.TW.func.do.set_func_body_as(velocity)
    # f1o.TW.current_time = 0
    # f1o.TW.___DO_push_all_to_instant___()
    # f1o.discretize()
    # f2 = FC('2-f-o', is_hybrid=False, name='total pressure')
    # f2.TW.func.do.set_func_body_as(scalar)
    # f2.TW.current_time = 0
    # f2.TW.___DO_push_all_to_instant___()
    # f2.discretize()
    #
    # grid = [np.linspace(-1,1,50), np.linspace(-1,1,50)]
    #
    # unstructuredGridToVTK(grid, [f0, f1o, f1i, f2], 'unstructuredGridToVTK_2dtest')