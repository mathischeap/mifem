# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, the Netherlands

"""
from abc import ABC
from time import localtime, strftime

from screws.miscellaneous.checks.check_no_splcharacter import check_no_splcharacter
from screws.miscellaneous.checks.check_multiple_close import check_multiple_close
from screws.miscellaneous.checks.check_filename import check_filename
from screws.miscellaneous.checks.check_filename_mi import check_filename_mi
from screws.miscellaneous.checks.check_almost_in_range import check_almost_in_range

from screws.miscellaneous.docstringReaders.numpy_styple import NumpyStyleDocstringReader
from screws.miscellaneous.randomString.digits import randomStringDigits
from screws.miscellaneous.break_list import break_list_into_parts
from screws.miscellaneous.get_my_IP_data import get_my_IP_data
from screws.miscellaneous.initialize_3d_list import initialize_3d_list
from screws.miscellaneous.count_files_and_lines import count_files_and_lines


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





if __name__ == '__main__':
    import doctest
    doctest.testmod()

    count_files_and_lines('../../')

    a = check_no_splcharacter
    b = check_multiple_close
    c = check_filename
    d = check_filename_mi
    e = check_almost_in_range
    f = NumpyStyleDocstringReader
    g = randomStringDigits
    h = break_list_into_parts
    i = initialize_3d_list
    IP = get_my_IP_data
