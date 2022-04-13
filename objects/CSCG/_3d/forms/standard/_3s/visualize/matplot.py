

from objects.CSCG._3d.forms.standard.base.visualize.matplot import _3dCSCG_standard_form_Matplot


class _3dCSCG_S3F_VISUALIZE_Matplot(_3dCSCG_standard_form_Matplot):
    """"""
    def __init__(self, sf):
        """"""
        super(_3dCSCG_S3F_VISUALIZE_Matplot, self).__init__(sf)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        raise NotImplementedError()