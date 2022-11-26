# -*- coding: utf-8 -*-

from components.miscellaneous.checks.check_no_splcharacter import check_no_splcharacter


def check_filename(filename):
    """A filename (with or without extension, if with, do not care about the exact extension) checker for files.

    .. doctest::

        >>> check_filename('a....12d')
        ('a___', '12d')
        >>> check_filename('__aba_dSSF.mii')
        ('__aba_dSSF', 'mii')

    """
    assert isinstance(filename, str) and check_no_splcharacter(filename), \
        f"filename={filename} is a {filename.__class__.__name__} wrong, need be a str."
    assert len(filename) <= 64, f"filename={filename} too long."

    if filename.count('.') > 1:
        str_list = filename.split('.')
        name_list = str_list[:-1]
        name = '_'.join(name_list)
        ext = str_list[-1]
        filename = name + '.' + ext

    assert filename.count('.') <= 1, f"filename={filename} is wrong: Can not have more than one 'dot' in the filename"

    if filename.count('.') == 0:
        file_name = filename
        extension = None
    elif filename.count('.') == 1:
        file_name, extension = filename.split('.')
    else:
        raise Exception(f'filename={filename} is illegal.')

    DIGITS = '0123456789'
    if extension is not None:
        ALL_digit = list()
        for e in extension:
            if e in DIGITS:
                ALL_digit.append(True)
            else:
                ALL_digit.append(False)

        if all(ALL_digit):
            raise Exception(f'extension={extension} is illegal, cannot be all digit.')

    HAVE_letter = list()
    LETTERS = 'qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM0123456789'
    for f in file_name:
        if f in LETTERS:
            HAVE_letter.append(True)
            break
        else:
            HAVE_letter.append(False)

    assert any(HAVE_letter), f"file_name={file_name} is illegal, need have a letter."
    assert file_name[0] not in DIGITS, f"file_name={file_name} is illegal, cannot start with a digit."
    assert ' ' not in file_name, f"file_name={file_name} contains space, wrong."

    return file_name, extension