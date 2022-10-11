# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/27/2022 12:30 AM
"""


from objects.miUsGrid.triangular.mesh.main import miUsGrid_TriangularMesh as mesh

from objects.miUsGrid.triangular.master import Call as form

call = form # another name of FormCaller

from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace as space