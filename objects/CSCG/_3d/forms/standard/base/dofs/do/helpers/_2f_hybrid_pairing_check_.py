# -*- coding: utf-8 -*-

from root.config.main import RANK, MASTER_RANK


def _2f_hybrid_pairing_check_(sf, adt, T, D, b):
    """

    Parameters
    ----------
    sf
    adt
    T
    D
    b

    Returns
    -------

    """
    assert str(adt.__class__.__name__) == "_3dCSCG_T2_ADF"
    assert sf.mesh == adt.mesh
    assert sf.space == adt.space

    Ta = T.assembled
    Da = D.assembled
    ba = b.assembled

    Ta = Ta.do.gather_M_to_core()
    Da = Da.do.gather_M_to_core()
    ba = ba.do.gather_V_to_core()

    tGB_dofs = adt.prime.numbering.GLOBAL_boundary_dofs
    sGB_dofs = sf.numbering.GLOBAL_boundary_dofs

    if RANK == MASTER_RANK:
        Ta = Ta.tocsr() # just in case
        Da = Da.tocsr() # just in case

        DB = adt.BC.boundaries
        NB = sf.BC.boundaries

        DB_dofs = set()
        for db in DB:
            DB_dofs.update(tGB_dofs[db])

        NB_dofs = set()
        for nb in NB:
            NB_dofs.update(sGB_dofs[nb])

        for i, Ti in enumerate(Ta):
            Ti_nnz = Ti.nnz
            if Ti_nnz == 2: # internal pairing
                assert Da[i].nnz == 0 and ba[i] == 0, f"Must be an internal dof (pair)."
                assert sum(Ti.data) == 0, f"it must be a pairing."
            elif Ti_nnz == 1: # on Neumann boundary
                data = Ti.data[0]
                assert data == 1, f"to impose Neumann boundary condition."
                assert ba[i] != 0, f"must have value in b."
                ind = Ti.indices[0]
                assert ind in NB_dofs, f"must be in the Neumann dofs."
                NB_dofs.remove(ind)
            elif Ti_nnz == 0: # on Dirichlet boundary
                assert Da[i].nnz == 1 and Da[i,i] == 1, f"to impose Dirichlet boundary condition."
                assert ba[i] != 0, f"must have value in b."
                assert i in DB_dofs, f"must be in the Dirichlet dofs."
                DB_dofs.remove(i)
            else:
                raise Exception()

        assert len(DB_dofs) == 0, f"We must have found (then removed) all Dirichlet dofs."
        assert len(NB_dofs) == 0, f"We must have found (then removed) all Neumann dofs."
