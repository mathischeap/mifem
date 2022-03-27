

from root.read.helpers._3dCSCG.form import ___restore__3dCSCG_Form___
from objects.CSCG._3d.master import FormCaller as _3dCSCG_FormCaller

def ___restore__3dCSCG_Algebra_DUAL_Form___(parameters, mesh_cache, space_cache):
    """
    For dual forms, we actually directly use the parameters from its prime,
    and the form type is decided by the obj_name.

    """
    assert parameters['type'] == '_3dCSCG_Form' # Not an error as we actually restore its prime form.
    ID = parameters['ID']
    prime = ___restore__3dCSCG_Form___(parameters, mesh_cache, space_cache)
    mesh = prime.mesh
    space = prime.space

    if ID == '0-f':
        dID = '0-adf'
    elif ID == '1-f':
        dID = '1-adf'
    elif ID == '2-f':
        dID = '2-adf'
    elif ID == '3-f':
        dID = '3-adf'
    elif ID == '0-t':
        dID = '0-adt'
    elif ID == '1-t':
        dID = '1-adt'
    elif ID == '2-t':
        dID = '2-adt'
    else:
        raise Exception()

    _3FC = _3dCSCG_FormCaller(mesh, space)

    dual = _3FC(dID, prime)

    return dual