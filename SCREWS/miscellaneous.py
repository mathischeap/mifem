# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from abc import ABC
from types import FunctionType, MethodType
import os
from time import localtime, strftime
import random
import string
from SCREWS.decorators import accepts
import re
import json
from urllib.request import urlopen

class NumpyStyleDocstringReader(ABC):
    """ This class can only read numpy style docstring."""
    def __init__(self, fm2read):
        isinstance(fm2read, (FunctionType, MethodType)), " <DocstringReader> : I need a function or method to read."
        self._fm2read_ = fm2read
        self._docstring_ = fm2read.__doc__
        self._docstring_no_space_ = self._docstring_.replace(' ', '')
        self.___Parameters___ = None
        self.___Returns___ = None
        self.___AllowedLinearGlobalSystem___ = None

    @property
    def Parameters(self):
        if self.___Parameters___ is None:
            _Para_ = self._docstring_no_space_.split('Parameters')[1]

            if '：' in _Para_:
                _Para_ =_Para_.replace('：', ":")

            _Para_ = _Para_.split('----------')[1]
            _Para_ = _Para_.split('\n\n')[0]
            _paras_ = _Para_.split(':')
            para_names = ()
            for ps in range(1, len(_paras_)):
                para_names += (_paras_[ps - 1].split('\n')[-1],)
            self.___Parameters___ = para_names
        return self.___Parameters___

    @property
    def Returns(self):
        if self.___Returns___ is None:
            _Re_ = self._docstring_no_space_.split('Returns')[1]

            if '：' in _Re_:
                _Re_ = _Re_.replace('：', ":")

            _Re_ = _Re_.split('-------')[1]
            _Re_ = _Re_.split('\n\n')[0]
            _Res_ = _Re_.split(':')
            Re_names = ()
            for rs in range(1, len(_Res_)):
                Re_names += (_Res_[rs - 1].split('\n')[-1],)
            self.___Returns___ = Re_names
        return self.___Returns___



class MyTimer(ABC):
    """My timers."""
    @classmethod
    def current_time(cls):
        """(str) Return a string showing current time."""
        return strftime("[%Y-%m-%d %H:%M:%S]", localtime())

    @classmethod
    def current_time_with_no_special_characters(cls):
        ct = cls.current_time()
        ct = ct.replace(' ', '_')
        ct = ct.replace('[', '')
        ct = ct.replace(']', '')
        ct = ct.replace(':', '_')
        ct = ct.replace('-', '_')
        return ct

    @classmethod
    def simple_header(cls):
        print("\n\t\t _____" + "___________________" + "_____")
        print("\t\t |^_^ " + " <TUD>-<YZ>-<AERO> " + " ^_^|")
        print("\t\t |>>> " + strftime("%Y-%m-%d %H:%M:%S", localtime()) + " <<<|\n")

    @classmethod
    def seconds2hms(cls, seconds):
        """We convert float: seconds to str: '[hh:mm:ss]'."""
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        return '[%02d:%02d:%02d]' % (hours, minutes, seconds)

    @classmethod
    def seconds2hmsm(cls, seconds):
        """We convert float: seconds to str: '[hh:mm:ss.ms]'."""
        ms = (seconds - int(seconds)) * 1000
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        return '[%02d:%02d:%02d.%03d]' % (hours, minutes, seconds, ms)

    @classmethod
    def seconds2dhms(cls, seconds):
        """We convert float: seconds to str: '[dd:hh:mm:ss]'."""
        days = int(seconds / 86400)
        SECONDS = seconds - 86400 * days
        minutes, SECONDS = divmod(SECONDS, 60)
        hours, minutes = divmod(minutes, 60)
        return '[{}:'.format(days) + '%02d:%02d:%02d]' % (hours, minutes, SECONDS)

    @classmethod
    def seconds2dhmsm(cls, seconds):
        """We convert float: seconds to str: '[dd:hh:mm:ss]'."""
        ms = (seconds - int(seconds)) * 1000
        days = int(seconds / 86400)
        SECONDS = seconds - 86400 * days
        minutes, SECONDS = divmod(SECONDS, 60)
        hours, minutes = divmod(minutes, 60)
        return '[{}:'.format(days) + '%02d:%02d:%02d.%03d]' % (hours, minutes, SECONDS, ms)

    @classmethod
    def hms2seconds(cls, hms):
        """We convert str: '[hh:mm:ss]' into float: seconds."""
        hh, mm, ss = hms[1:-1].split(':')
        hh = int(hh) * 3600
        mm = int(mm) * 60
        ss = int(ss)
        return hh + mm + ss

    @classmethod
    def dhms2seconds(cls, hms):
        """We convert str: '[hh:mm:ss]' into float: seconds."""
        dd, hh, mm, ss = hms[1:-1].split(':')
        dd = int(dd) * 86400
        hh = int(hh) * 3600
        mm = int(mm) * 60
        ss = int(ss)
        return dd + hh + mm + ss

    @classmethod
    def dhmsm2seconds(cls, hmsm):
        """We convert str: '[hh:mm:ss]' into float: seconds."""
        hms, ms = hmsm.split('.')
        ms = int(ms[:-1]) / 1000
        dd, hh, mm, ss = hms[1:].split(':')
        dd = int(dd) * 86400
        hh = int(hh) * 3600
        mm = int(mm) * 60
        ss = int(ss)
        return dd + hh + mm + ss + ms



def randomStringDigits(stringLength=8):
    """Generate a random string of letters and digits."""
    lettersAndDigits = string.ascii_letters + string.digits
    return ''.join(random.choice(lettersAndDigits) for _ in range(stringLength))



def initialize_3d_list(a, b, c):
    """

    :param a:
    :param b:
    :param c:
    :return:
    """
    lst = [[[None for _ in range(c)] for _ in range(b)] for _ in range(a)]
    return lst


@accepts((int, float), (int, float))
def check_multiple_close(a, b, tol=1e-8):
    """check if a = b*i +- tol where i = 1,2,3,4,...

    :param a:
    :param b:
    :param tol:
    :return:
    """
    remainder = a % b
    if remainder < tol:
        return True
    else:
        assert b > remainder, "something wrong."
        if (b - remainder) < tol:
            return True
    return False


@accepts((int, float), (int, float), (int, float))
def check_almost_in_range(a, lb, ub, tol=1e-8):
    """check a is in [lb-tol, up+tol].

    :param a:
    :param lb:
    :param ub:
    :param tol:
    :return:
    """
    assert lb < ub, f"lower bound {lb} must be lower than upper bound {ub}."
    return (lb - tol) < a < (ub + tol)


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


@accepts(str)
def check_no_splcharacter(test):
    """Check the string ``test`` have no special character listed below.

    :param test:
    :return:
    """
    # noinspection RegExpRedundantEscape
    string_check = re.compile('[@!#$%^&*()<>?/\|}{~:]')
    if string_check.search(test) is None:
        return True
    else:
        return False


def get_my_IP_data():
    """

    :return:
    """
    url = 'http://ipinfo.io/json'
    response = urlopen(url)
    data = json.load(response)
    return data


def break_list_into_parts(LIST, parts):
    """For example:

    .. doctest::

        >>> break_list_into_parts([1,2,3,4,5], [3,2])
        ([1, 2, 3], [4, 5])

    :param LIST:
    :param parts:
    :return:
    """
    start = 0
    RET = tuple()
    for p in parts:
        RET += (LIST[start:start+p],)
        start += p
    return RET


def count_files_and_lines(start, files=0, lines=0, header=True, begin_start=None):
    """
    We use this function to count how many python files are in a path and how
    many lines of python code are in these files.

    The original code is provided by stack overflow user 'Bryce93' on page
    'https://stackoverflow.com/questions/38543709/count-lines-of-code-in-directory-using-python'.

    We have skipped
        1). files in folder '__PROGRAMS__'
        2). meaningless lines, i.e. of only spaces or '-'.
        3). files named with extension '.pyc'.

    Parameters
    ----------
    start : str
        The path we want to start tracing in.
    files :
    lines :
    header :
    begin_start :

    """
    if header:
        MyTimer.simple_header()
        print("\t\t Line counting......\n")
        print('{:>10} |{:>10} | {:<20}'.format('ADDED', 'TOTAL', 'FILE'))
        print('{:->11}|{:->11}|{:->20}'.format('', '', ''))
    toskip = ['    ' * i + '\n' for i in range(12)]
    for thing in os.listdir(start):
        thing = os.path.join(start, thing)
        if os.path.isfile(thing):
            if thing.endswith('.py'): # we only look at python files.
                with open(thing, 'r') as f:
                    try:
                        FILES = f.readlines()[9:]
                    except UnicodeDecodeError:
                        lines += 20
                    else:
                        newlines = len(FILES)
                        for entry in FILES:
                            if entry in toskip:
                                newlines -= 1
                            elif '____' in entry or '-' in entry:
                                newlines -= 1
                            elif '~' in entry or '===' in entry:
                                newlines -= 1
                            elif '"""' in entry or '%%' in entry:
                                newlines -= 1
                            else:
                                pass
                        newlines = int(newlines * 0.95) + 1  # we consider 5% of lines are docstring.
                        lines += newlines
                        if begin_start is not None:
                            reldir_of_thing = '.' + thing.replace(begin_start, '')
                        else:
                            reldir_of_thing = '.' + thing.replace(start, '')
                        if len(reldir_of_thing) > 45:
                            reldir_of_thing = reldir_of_thing[:42] + '...'
                        print('{:>10} |{:>10} | {:<20}'.format(newlines, lines, reldir_of_thing))
                        files += 1

    for thing in os.listdir(start):
        # skip following ...
        if thing[-4:] == '.pyc':
            pass
        elif thing[:12] == '_CONTENTS_':
            pass
        elif thing[:11] == '__GENERAL__':
            pass
        elif thing[:4] == 'venv':
            pass
        else:
            thing = os.path.join(start, thing)
            if os.path.isdir(thing):
                files, lines = count_files_and_lines(thing, files, lines, header=False, begin_start=start)
    return files, lines



if __name__ == '__main__':
    # a = check_multiple_close(1.9999999999999, 1)
    # print(a)
    # a = check_almost_in_range(2.000001, 1,2)
    # print(a)
    # # Enter the string to be checked
    #
    # test = "Code_Speedya adsas afeas @@@gasd"
    #
    # # calling check_splcharacter function
    #
    # a = check_no_splcharacter(test)
    # print(a)

    # SendAdminAnEmail('TEST message')()
    # from config import usErs
    # SendAdminAnHTMLEmail(usErs)('MESSAGE test')

    import doctest
    doctest.testmod()

    count_files_and_lines('./')