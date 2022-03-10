

from screws.miscellaneous.checks.check_no_splcharacter import check_no_splcharacter

def check_filename_mi(filename):
    """A special filename checker for '.mi' files."""
    assert isinstance(filename, str) and check_no_splcharacter(filename), \
        f"filename={filename} is a {filename.__class__.__name__} wrong, need be a str."
    assert len(filename) <= 64, f"filename={filename} too long."
    if filename[-3:] != '.mi': filename += '.mi'
    if filename.count('.') > 1:
        filename = filename.replace('.','_')
        filename = filename[:-3] + '.mi'
    assert filename.count('.') == 1 and filename[-3:] == '.mi', f"filename={filename} is wrong."
    assert ' ' not in filename, f"filename={filename} contains space, wrong."
    return filename