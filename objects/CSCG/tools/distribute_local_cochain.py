# -*- coding: utf-8 -*-
"""If we store local cochains of multiple forms in a dict, we can use this function to distribute
the local cochain to each form.

"""

def distribute_local_cochain(local_cochain, forms):
    """

    Parameters
    ----------
    local_cochain : dict
        A dict whose keys are local mesh elements, and values are the concatenated local cochains
        of `forms`.

    forms : list, tuple
        The forms who are going to receive the local cochains.

    Returns
    -------

    """
    start = 0

    for f in forms:

        num_basis = f.num.basis

        end = start + num_basis

        f_l_c = dict() # form local cochain

        for me in local_cochain:
            f_l_c[me] = local_cochain[me][start : end]

        f.cochain.local = f_l_c

        start += num_basis