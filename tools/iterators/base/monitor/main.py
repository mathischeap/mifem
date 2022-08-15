# -*- coding: utf-8 -*-
from root.config.main import *
from screws.freeze.main import FrozenOnly
import socket
import datetime
from screws.miscellaneous.timer import MyTimer
from time import time
from tools.iterators.base.monitor.do import IteratorMonitorDo
from tools.iterators.base.monitor.IS import IteratorMonitorIS

class IteratorMonitor(FrozenOnly):
    """The monitor class for the Iterator.

    :param iterator:
    :param auto_save_frequency:
    :param RDF_filename:
    :param factor:
    :param real_time_monitor:
    """
    def __init__(self, iterator, auto_save_frequency, RDF_filename, factor: float, real_time_monitor: bool):
        assert rAnk == mAster_rank, "Should only initialize and use it in master core."
        self.___graph_report_default_time___ = 600 # seconds; 10 minutes     --> affected by ``factor``
        self.___auto_save_default_time___ = 600  # seconds: 10 minutes       --> affected by ``factor``
        self.___email_report_default_time___ = 3600 * 12 # seconds; 12 hours --> affected by ``factor``
        self.___PRIVATE_set_factor___(factor)
        self.___PRIVATE_set_auto_save_frequency___(auto_save_frequency)
        self.___PRIVATE_set_RDF_filename___(RDF_filename)
        self._real_time_monitor_ = real_time_monitor

        self._iterator_ = iterator
        # Following are all timing related ...
        if self._iterator_._max_steps_ is not None and self._iterator_._max_steps_ > 999999:
            self._isOpen_ = True # the iterator will stop when recv shut_down order.
        else:
            self._isOpen_ = False
        self._max_steps_ = self._iterator_.max_steps
        self._times_ = np.array([np.nan,]) # cost of all iterations
        self._TIMES_ = np.array([0.,]) # total cost after each iteration.

        # the time initialize the iterator, mostly immediately after starting the program.
        self._ft_start_time_ = time()
        self._ft_firstRun_ = None # to be updated once when run starts.
        self._str_started_time_ = MyTimer.current_time()[1:-1]
        self._do_first_assessment_ = True
        self._do_first_assessment_counter_ = 0
        self._do_first_assessment_estimate_ = 0
        self._do_first_email_warning_report_ = False
        self._do_first_graph_warning_report_ = False
        self._ever_do_email_report_ = False
        self._ever_do_graph_report_ = False
        # to be updated after every iteration.
        self._computed_steps_ = 0
        self._preparation_time_ = None
        self._current_time_ = None
        self._total_cost_ = 0
        self._last_run_cost_ = None
        self._effective_run_cost_ = 0
        self._effective_run_num_ = 0
        self._average_each_run_cost_ = 1
        self._estimated_remaining_time_ = 1
        self._estimated_end_time_ = datetime.datetime.now()
        # ...
        self._do_ = IteratorMonitorDo(self)
        self._IS_ = IteratorMonitorIS(self)
        self._freeze_self_()

    def ___PRIVATE_set_factor___(self, factor):
        """"""
        if factor in (False, 0):
            factor = 0
        elif factor is True:
            factor = 1
        else:
            assert factor >= 0.1, \
                f"monitor factor={factor} wrong, should be `False`, `True` or `>= 0.1`"
        if factor == 0:
            self.___graph_report_time___ = 9999999
            self.___auto_save_time___    = 9999999
            self.___email_report_time___ = 9999999
        else:
            self.___graph_report_time___ = (1 / factor) * self.___graph_report_default_time___
            self.___auto_save_time___ = (1 / factor) * self.___auto_save_default_time___
            self.___email_report_time___ = (1 / factor) * self.___email_report_default_time___

        if self.___graph_report_time___ < 60: self.___graph_report_time___ = 60   # 最快60秒报告一次
        if self.___auto_save_time___ < 60: self.___auto_save_time___ = 60         # 最快60秒报告一次
        if self.___email_report_time___ < 3600: self.___email_report_time___ = 3600 # 最快1小时报告一次

        self._factor_ = factor
        # even it is 0, the monitor still do the recording background, but no report.
        self.___last_graph_save_time___ = time()
        self.___last_email_sent_time___ = time()


    def ___PRIVATE_set_auto_save_frequency___(self, AS):
        self._last_auto_save_time_ = None
        if AS is True:
            # will parse ``DO_auto_save`` after first iteration.
            self._auto_save_frequency_ = True
        elif AS is False:
            self._auto_save_frequency_ = 0
        elif isinstance(AS, (int, float)):
            assert AS >= 0, \
                " <iterator> : DO_auto_save frequency={} needs >= 0!".format(AS)
            self._auto_save_frequency_ = int(AS)
        else:
            raise Exception(" <iterator> : DO_auto_save frequency wrong!")
    def ___PRIVATE_set_RDF_filename___(self, RDF_filename):
        if RDF_filename is None:
            pass
        elif isinstance(RDF_filename, str):
            if '.' not in RDF_filename: RDF_filename += '.csv'
            assert RDF_filename[-4:] == '.csv', f"Need to be a csv file!, now it is {RDF_filename}."
        else:
            raise Exception("RDF_filename must be str. Now, it is %r."%RDF_filename)
        self._RDF_filename_ = RDF_filename

    @property
    def factor(self):
        """
        (float) Monitor factor, be in ``[0,1]``. When it is ``0``, no monitor. When it is ``1``,
        very frequent monitor.
        """
        return self._factor_

    @property
    def auto_save_frequency(self):
        """(int) The frequency of doing the auto-saving. """
        return self._auto_save_frequency_

    @property
    def RDF_filename(self):
        """(str)"""
        return self._RDF_filename_

    @property
    def IS(self):
        return self._IS_

    @property
    def summary_html(self):
        local_IP = socket.gethostbyname(socket.gethostname())
        local_machine_name = socket.gethostname()
        if self._estimated_remaining_time_ <= 0.1:
            EET = 'NOW'
        else:
            EET = str(self._estimated_end_time_)[:19]
        summary_html = \
        f"""
        <h2 style="background-color:PaleTurquoise;">
        {self._iterator_.__class__.__name__}: [{
            self._iterator_.standard_properties.name}] is running on [{local_machine_name}] @ [{local_IP}].
        </h2>

        <h2>Summery:</h2>

        <table>
          <tr>
            <th>item</th>
            <th>message</th>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-completed/total steps </b></em></td>
            <td>{self._computed_steps_}/{self._max_steps_}</td>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-started at </b></em></td>
            <td><font color="Crimson">{self._str_started_time_}<font></td>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-preparation cost </b></em></td>
            <td><font color="DarkBlue">{MyTimer.seconds2dhmsm(self._preparation_time_)}<font></td>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-total cost </b></em></td>
            <td><font color="DarkBlue">{MyTimer.seconds2dhmsm(self._total_cost_)}<font></td>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-iterations cost </b></em></td>
            <td>{MyTimer.seconds2dhmsm(self._effective_run_cost_)}</td>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-last iteration cost </b></em></td>
            <td>{MyTimer.seconds2dhmsm(self._last_run_cost_)}</td>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-effective iteration average cost </b></em></td>
            <td>{MyTimer.seconds2dhmsm(self._average_each_run_cost_)}</td>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-estimated remaining time</b></em></td>
            <td><font color="DarkBlue">{MyTimer.seconds2dhmsm(self._estimated_remaining_time_)}<font></td>
          </tr>
          <tr>
            <td style="border:2px solid DodgerBlue;"><em><b>-estimated completion at</b></em></td>
            <td><font color="Crimson">{EET}</font></td>
          </tr>
        </table>

        <h2>Solver message:</h2> 

        """
        for message in self._iterator_.message:
            summary_html += f"""
            <p>{message}</p>

            """
        summary_html += """
        <h2>Results:</h2> 

        """
        return summary_html

    @property
    def do(self):
        return self._do_
