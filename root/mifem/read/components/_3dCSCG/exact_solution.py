

from _3dCSCG.master import ExactSolutionSelector as _3dCSCG_ExactSolutionSelector
from root.mifem.read.components._3dCSCG.mesh import ___restore__3dCSCG_Mesh___





def ___restore__3dCSCG_ExactSolution___(parameters, mesh_cache):
    assert parameters.pop('type') == '_3dCSCG_ExactSolution'
    ID = parameters.pop('ID')
    mesh_parameters = parameters.pop('mesh_parameters')
    kwargs = parameters.pop('kwargs')
    assert len(parameters) == 0, "make sure all information are used."
    mesh = ___restore__3dCSCG_Mesh___(mesh_parameters, mesh_cache)
    return _3dCSCG_ExactSolutionSelector(mesh)(ID, **kwargs)