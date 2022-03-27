

from objects.CSCG._2d.master import ExactSolutionSelector as _2dCSCG_ExactSolutionSelector
from root.read.helpers._2dCSCG.mesh import ___restore__2dCSCG_Mesh___


def ___restore__2dCSCG_ExactSolution___(parameters, mesh_cache):
    assert parameters.pop('type') == '_2dCSCG_ExactSolution'
    ID = parameters.pop('ID')
    mesh_parameters = parameters.pop('mesh_parameters')
    kwargs = parameters.pop('kwargs')
    assert len(parameters) == 0, "make sure all information are used."
    mesh = ___restore__2dCSCG_Mesh___(mesh_parameters, mesh_cache)
    return _2dCSCG_ExactSolutionSelector(mesh)(ID, **kwargs)
