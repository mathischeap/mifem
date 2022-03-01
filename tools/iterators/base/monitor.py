
import os
from root.config import *
from screws.frozen import FrozenOnly
import socket
import codecs
from screws.emails import SendAdminAnHTMLEmail, whether_internet_connected, SendAdminAnEmail
import matplotlib.pyplot as plt
from matplotlib import cm
import datetime
from screws.miscellaneous import MyTimer
from time import time, sleep



class IteratorMonitor(FrozenOnly):
    """The monitor class for the Iterator.

    :param iterator:
    :param auto_save_frequency:
    :param RDF_filename:
    :param factor:
    """
    def __init__(self, iterator, auto_save_frequency, RDF_filename, factor: float):
        assert rAnk == mAster_rank, "Should only initialize and use it in master core."
        self.___auto_save_default_time___ = 3600 # seconds: 1 hour           --> affected by ``factor``
        self.___graph_report_default_time___ = 600 # seconds; 10 minutes     --> affected by ``factor``
        self.___email_report_default_time___ = 3600 * 12 # seconds; 12 hours --> affected by ``factor``
        self.___PRIVATE_set_factor___(factor)
        self.___PRIVATE_set_auto_save_frequency___(auto_save_frequency)
        self.___PRIVATE_set_RDF_filename___(RDF_filename)

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
        self._freeze_self_()

    @property
    def factor(self):
        """
        (float) Monitor factor, be in ``[0,1]``. When it is ``0``, no monitor. When it is ``1``,
        very frequent monitor.
        """
        return self._factor_
    def ___PRIVATE_set_factor___(self, factor):
        if factor in (False, 0):
            factor = 0
        elif factor is True:
            factor = 1
        else:
            assert factor >= 0.1, f"monitor factor={factor} wrong, should be False, True or >=0.1"
        if factor == 0:
            self.___graph_report_time___ = 999999
            self.___email_report_time___ = 999999
            self.___auto_save_time___    = 999999
        else:
            self.___graph_report_time___ = (1 / factor) * self.___graph_report_default_time___
            self.___email_report_time___ = (1 / factor) * self.___email_report_default_time___
            self.___auto_save_time___ = (1 / factor) * self.___auto_save_default_time___

        if self.___graph_report_time___ < 60: self.___graph_report_time___ = 60
        if self.___email_report_time___ < 720: self.___email_report_time___ = 720
        if self.___auto_save_time___ < 600: self.___auto_save_time___ = 600

        self._factor_ = factor
        # even it is 0, the monitor still do the recording background, but no report.
        self.___last_graph_save_time___ = time()
        self.___last_email_sent_time___ = time()

    @property
    def auto_save_frequency(self):
        """(int) The frequency of doing the auto-saving. """
        return self._auto_save_frequency_
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

    @property
    def RDF_filename(self):
        """(str)"""
        return self._RDF_filename_
    def ___PRIVATE_set_RDF_filename___(self, RDF_filename):
        if RDF_filename is None:
            pass
        elif isinstance(RDF_filename, str):
            if '.' not in RDF_filename: RDF_filename += '.csv'
            assert RDF_filename[-4:] == '.csv', f"Need to be a csv file!, now it is {RDF_filename}."
        else:
            raise Exception("RDF_filename must be str. Now, it is %r."%RDF_filename)
        self._RDF_filename_ = RDF_filename

    def ___PRIVATE_select_reasonable_amount_of_data___(self, max_num, last_num=1):
        """
        To report RDF, we do not report all, we make a selection.

        :param max_num: we at most selection this amount of data
        :param last_num: we must select the last ``last_num`` rows of RDF.
        """
        assert max_num >= 10
        assert 1 <= last_num < max_num
        all_data_num = len(self._iterator_.RDF)
        if max_num < all_data_num:
            indices = list(np.linspace(0, all_data_num-last_num, max_num+1-last_num).astype(int))
            indices.extend([all_data_num+1-last_num+i for i in range(last_num-1)])
        else:
            indices = [i for i in range(all_data_num)]
        return indices

    @property
    def IS_open(self):
        """(bool) Return ``True`` if the iterator is open (no max_steps, stop when shut_down)."""
        return self._isOpen_

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

    def DO_update(self):
        """Update monitor after every iteration."""
        self._computed_steps_ += 1
        if self._current_time_ is None:
            self._current_time_ = self._ft_firstRun_ # to make the first iteration time correct.
        self._last_run_cost_ = time() - self._current_time_
        self._current_time_ = time()
        self._total_cost_ = self._current_time_ - self._ft_start_time_
        if self._last_run_cost_ > 0.25:  # we only consider cost long enough iteration as effective.
            self._effective_run_cost_ += self._last_run_cost_
            self._effective_run_num_ += 1
            self._average_each_run_cost_ = self._effective_run_cost_ / self._effective_run_num_
            self._estimated_remaining_time_ = self._average_each_run_cost_ * (
                    self._max_steps_ - self._computed_steps_
            )
            if self._estimated_remaining_time_ > 3600*24*99:
                self._estimated_remaining_time_ = 3600*24*99 + 23*3600 + 59*60 + 59.999
            self._estimated_end_time_ = datetime.datetime.now() + datetime.timedelta(
                seconds=self._estimated_remaining_time_)

            if self._do_first_assessment_: # will see first 5 iterations to see the estimated time.
                self._do_first_assessment_counter_ += 1
                self._do_first_assessment_estimate_ += self._estimated_remaining_time_
                est = self._do_first_assessment_estimate_ / self._do_first_assessment_counter_
                if est > self.___email_report_time___:
                    self._do_first_email_warning_report_ = True
                if est > self.___graph_report_time___:
                    self._do_first_graph_warning_report_ = True
                if (self._do_first_email_warning_report_ and self._do_first_graph_warning_report_) or \
                    self._do_first_assessment_counter_ >= 5:
                    self._do_first_assessment_ = False
                    del self._do_first_assessment_counter_, self._do_first_assessment_estimate_

            self._times_ = np.append(self._times_, self._last_run_cost_)
        else:
            # this is a interesting choice. No clue why I chose do so
            self._times_ = np.append(self._times_, np.nan)
        self._TIMES_ = np.append(self._TIMES_, self._total_cost_)

        if (self._max_steps_ == self._computed_steps_) or self._iterator_.shut_down:
            # make sure we have correct _estimated_remaining_time_ when iterator is done
            self._estimated_remaining_time_ = 0
            self._estimated_end_time_ = datetime.datetime.now()

    # noinspection PyBroadException
    def DO_auto_save(self):
        if self.RDF_filename is not None:
            if self._last_auto_save_time_ is None:
                self._last_auto_save_time_ = self._ft_firstRun_
            gap_time = time() - self._last_auto_save_time_
            if self.auto_save_frequency is True:
                if gap_time > self.___auto_save_time___:
                    try:
                        # if PermissionError, we do not stop the iteration
                        self._iterator_.RDF.to_csv(self.RDF_filename, header=True)
                        self._last_auto_save_time_ = time()
                    except: # wait 10 seconds
                        sleep(10)
                        try: # try once more
                            self._iterator_.RDF.to_csv(self.RDF_filename, header=True)
                        except: # just skip it
                            pass
                else:
                    pass
            elif self.auto_save_frequency > 0:
                if ((self._computed_steps_ % self.auto_save_frequency == 0) and
                    gap_time > 0.1*self.___auto_save_time___) or \
                    gap_time > self.___auto_save_time___:
                    # important. When read from csv, it is very fast, so we do not save.
                    try: # if PermissionError, we do not stop the iteration
                        self._iterator_.RDF.to_csv(self.RDF_filename, header=True)
                        self._last_auto_save_time_ = time()
                    except: # wait 10 seconds
                        sleep(10)
                        try: # try once moe
                            self._iterator_.RDF.to_csv(self.RDF_filename, header=True)
                        except: # just skip it
                            pass
            else:
                pass
        else:
            pass

    def ___PRIVATE_generate_open_graph_report___(self):
        """"""
        raise NotImplementedError()

    def ___PRIVATE_generate_open_email_report___(self):
        """"""
        raise NotImplementedError()

    def DO_generate_graph_report(self):
        """"""
        # do A LAST REPORT BEFORE STOP:
        judge1 = \
            ((self._computed_steps_ == self._max_steps_) or self._iterator_.shut_down) and \
            (self._total_cost_ > self.___graph_report_time___ or self._ever_do_graph_report_)

        # intermediate save: has cost a certain time or it will cost long time, so we do a first report
        judge2 = ((self._current_time_ - self.___last_graph_save_time___) > self.___graph_report_time___) or \
            self._do_first_graph_warning_report_

        # to judge special iteration stop: shut down by the solver.
        judge_sd = self._iterator_.shut_down

        if self._do_first_graph_warning_report_: self._do_first_graph_warning_report_ = False

        if judge1 or judge2: # now we need to do the reporting.

            if not self._ever_do_graph_report_: self._ever_do_graph_report_ = True

            if self.IS_open: return self.___PRIVATE_generate_open_graph_report___()

            save_time = MyTimer.current_time()[1:-1]
            indices = self.___PRIVATE_select_reasonable_amount_of_data___(1000, last_num=100)
            RDF = self._iterator_.RDF.iloc[indices]
            plt.rc('text', usetex=False)
            num_subplots = RDF.shape[1] + 3 # We plot 2 extra: 't iteration' and 't accumulation' + solver message
            colors = cm.get_cmap('cool_r', num_subplots - 5)
            r_num_subplots = int(np.ceil(num_subplots/2))
            x_len, y_len = 18, 4.5*r_num_subplots
            fig = plt.figure(figsize=(x_len, y_len))
            plot_keys = list(RDF.columns)
            plot_keys = plot_keys[:2] + ['t iteration', 't accumulation', 'solver message',] + plot_keys[2:]

            # subplots ...
            for i, di in enumerate(plot_keys):
                ylabel_backgroundcolor = 'paleturquoise'
                face_color = 'lavender'
                ylabel = di.replace('_', '-')
                m = int(i/2)
                n = i % 2
                # noinspection PyUnboundLocalVariable
                ax = plt.subplot2grid((r_num_subplots, 2),(m, n))
                if di == 'solver message':
                    face_color = "whitesmoke"
                    plt.axis([0, 10, 0, 10])
                    ax.axes.get_xaxis().set_visible(False)
                    ax.axes.get_yaxis().set_visible(False)
                    plt.text(0.1, 9.5, 'SOLVER MESSAGE:', color= 'deepskyblue', fontsize=18, style='normal',
                             ha='left', va='top', wrap=True)
                    message = ''
                    for M in self._iterator_.message:
                        if len(M) < 70:
                            message += M
                        else:
                            MS = M.split(' ')
                            new_line = 0
                            for ms in MS:
                                if ms == '':
                                    # remove extra space.
                                    pass
                                else:
                                    new_line += len(ms) + 1
                                    message += ms + ' '
                                    if new_line >= 60:
                                        message += '\n'
                                        new_line = 0
                        message += '\n\n'
                    plt.text(0.1, 8.5, message, color= 'black', fontsize=12,
                             ha='left', va='top', wrap=True)
                elif di == 't':
                    face_color = 'honeydew'
                    plt.axis([0, 10, 0, 10])
                    ax.axes.get_xaxis().set_visible(False)
                    ax.axes.get_yaxis().set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    ax.spines['left'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    t1 = f'* {self._iterator_.standard_properties.name}'
                    plt.text(-2, 11.3, t1, color= 'darkorchid', fontsize=26, style='normal', ha='left',
                             va='top', wrap=True)
                    sITC = MyTimer.seconds2hms(self._average_each_run_cost_)
                    sPPT = MyTimer.seconds2hms(self._preparation_time_)
                    sTTC = MyTimer.seconds2dhmsm(self._total_cost_) .split('.')[0]  + ']'
                    sLIC = MyTimer.seconds2hms(self._last_run_cost_)
                    sERT = MyTimer.seconds2dhmsm(self._estimated_remaining_time_).split('.')[0]  + ']'
                    t2 = 'ITC: ' + sITC  + '      LIC: '    + sLIC + '\n'
                    t2 += 'TTC: ' + sTTC + '  Preparation:' + sPPT + ' \n'
                    t2 += 'ERT: ' + sERT + '\n'
                    percentage = int(10000*(self._computed_steps_/self._max_steps_)) / 100
                    t2 += f'Iterations done: {self._computed_steps_}/{self._max_steps_} ~ {percentage}%\n'
                    t2 += f'Iterator type: {self._iterator_.__class__.__name__}\n'
                    plt.text(-0.5, 9.5, t2, color= 'darkblue', fontsize=22, style='normal', ha='left',
                             va='top', wrap=True)
                    t3 = 'Graph saved at:: ' + save_time + '\n'
                    t3 += 'Start at :: ' + self._str_started_time_ + '\n'
                    if self._estimated_remaining_time_ == 0:
                        t3 += 'E end at:: NOW'
                    else:
                        t3 += 'E end at:: ' + str(self._estimated_end_time_)[:19]
                    plt.text(-0.5, 3.25, t3, color= 'red', fontsize=22, style='normal', ha='left',
                             va='top', wrap=True)
                    local_IP = socket.gethostbyname(socket.gethostname())
                    local_machine_name = socket.gethostname()
                    t4 = "Running on <" + local_machine_name + '@' + local_IP +'>'
                    plt.text(-1, -0.5, t4, color= 'black', fontsize=20, style='normal', ha='left',
                             va='top', wrap=True)

                elif di == 'dt':
                    ylabel_backgroundcolor = 'greenyellow'
                    face_color = 'snow'

                    if len(indices) < 25:
                        plt.plot(indices, RDF[di], '-o', color='k', linewidth=0.8)
                    else:
                        plt.plot(indices, RDF[di], color='k', linewidth=1.2)
                    if len(indices) <= 5:
                        ax.set_xticks(indices)

                    ylim = ax.get_ylim()
                    plt.text(0, ylim[1] - (ylim[1] - ylim[0])*0.1, "saved at: "+save_time, color= 'teal',
                             fontsize=22, style='normal', ha='left',
                             va='top', wrap=True)

                    from_t = self._iterator_.RDF['t'][indices[0]]
                    plt.text(0, ylim[0] + (ylim[1] - ylim[0])*0.2, f"compute from: t=%.3f s."%from_t, color= 'k',
                             fontsize=22, style='normal', ha='left',
                             va='bottom', wrap=True)
                    till_t = self._iterator_.RDF['t'][indices[-1]]
                    plt.text(0, ylim[0] + (ylim[1] - ylim[0])*0.1, f"compute till: t=%.3f s."%till_t, color= 'k',
                             fontsize=22, style='normal', ha='left',
                             va='bottom', wrap=True)

                elif di == 't iteration':
                    ylabel_backgroundcolor = 'greenyellow'
                    face_color = 'snow'
                    itime = self._times_[indices]
                    itime = self.___DO_filter_extreme_time___(itime) # extreme values are replaced by nan.
                    valid_time = itime[~np.isnan(itime)]
                    max_time = np.max(valid_time)
                    if max_time < 999:
                        ylabel = 't iteration (s)'
                        unit = 's'
                    elif max_time < 3600 * 3:
                        itime /= 60
                        ylabel = 't iteration (m)'
                        unit = 'm'
                    else:
                        itime /= 3600
                        ylabel = 't iteration (h)'
                        unit = 'h'
                    if len(indices) < 25:
                        plt.plot(indices, itime, '-o', color='k', linewidth=0.8, label='ITC')
                    else:
                        plt.plot(indices, itime, color='k', linewidth=1.2, label='ITC')

                    v10_ratio = None
                    if len(valid_time) >= 10 and np.ndim(valid_time) == 1:
                        average = np.sum(valid_time) / len(valid_time)

                        last_10_time = valid_time[-10:]
                        average_last_10 = np.sum(last_10_time) / 10

                        if isinstance(average_last_10, (int, float)) and average_last_10 > 0:
                            if isinstance(average, (int, float)) and average > 0:
                                v10_ratio = average_last_10 / average

                    else:
                        average = self._average_each_run_cost_
                        average_last_10 = None


                    if unit == 's':
                        pass
                    elif unit == 'm':
                        average /= 60
                        if average_last_10 is not None: average_last_10 /= 60
                    elif unit == 'h':
                        average /= 3600
                        if average_last_10 is not None: average_last_10 /= 3600
                    else:
                        raise Exception()
                    if indices[-1] > 1:
                        plt.plot([1, indices[-1]], [average, average],
                                 color='red', linewidth=1.2, label='AVR')

                    if average_last_10 is None:
                        pass
                    else:
                        plt.plot([1, indices[-10]], [average_last_10, average_last_10], '--',
                                 color='blue', linewidth=1.2)
                        plt.plot([indices[-10],indices[-1]], [average_last_10, average_last_10],
                                 color='blue', linewidth=1.2, label='V10')


                    if len(indices) <= 5: ax.set_xticks(indices)
                    plt.legend(fontsize=16, loc='best', frameon=False)

                elif di == 't accumulation':
                    ylabel_backgroundcolor = 'greenyellow'
                    face_color = 'snow'
                    TIME = self._TIMES_[indices]
                    MAX_TIME= TIME[-1]
                    if MAX_TIME < 60 * 3:
                        ylabel = 'Cumulative t (s)'
                    elif MAX_TIME < 3600 * 3:
                        TIME /= 60
                        ylabel = 'Cumulative t (m)'
                    elif MAX_TIME < 3600 * 24 * 3:
                        TIME /= 3600
                        ylabel = 'Cumulative t (h)'
                    else:
                        TIME /= 3600 * 24
                        ylabel = 'Cumulative t (d)'
                    if len(indices) < 25:
                        plt.plot(indices, TIME, '-o', color='k', linewidth=0.8)
                    else:
                        plt.plot(indices, TIME, color='k', linewidth=1.2)
                    if len(indices) <= 5:
                        ax.set_xticks(indices)

                    if not judge1:
                        sERT = MyTimer.seconds2dhmsm(self._estimated_remaining_time_).split('.')[0] + ']'
                        y_position = 0.9 * TIME[-1]
                        plt.text(0, y_position, 'ERT: '+ sERT, color= 'darkblue',
                                 fontsize=22, style='normal', ha='left',
                                 va='top', wrap=True)
                        # noinspection PyUnboundLocalVariable
                        if v10_ratio is not None:
                            # noinspection PyUnboundLocalVariable
                            v10_ERT_seconds = self._estimated_remaining_time_ * v10_ratio
                            vERT = MyTimer.seconds2dhmsm(v10_ERT_seconds).split('.')[0] + ']'
                            y_position = 0.75 * TIME[-1]
                            plt.text(0, y_position, 'V10: '+ vERT, color= 'purple',
                                     fontsize=22, style='normal', ha='left',
                                     va='top', wrap=True)

                else:
                    plt.plot(RDF['t'], RDF[di], color=colors(i-4), linewidth=1.5)

                ax.tick_params(axis="x", direction='in', length=8, labelsize=15)
                ax.tick_params(axis="y", direction='in', length=8, labelsize=15)
                # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0)) # set y ticks to sci format.
                tx = ax.yaxis.get_offset_text()
                # ... change the ytick sci-format font size
                tx.set_fontsize(15)
                ax.set_ylabel(ylabel, fontsize=17, backgroundcolor=ylabel_backgroundcolor)

                # ITERATING watermark ...
                if i < 4: # 0,1,2,3 for regular plots, no ITERATING
                    pass
                elif di == 'solver message': # no ITERATING for solver message subplot
                    pass
                elif not judge1: # if not the last iteration, texture it.
                    the_text = 'ITERATING'
                    text = ax.text(0.1, 0.5, the_text, fontsize=65, color='gray',
                                   horizontalalignment='left',
                                   verticalalignment='center',
                                   transform=ax.transAxes)
                    text.set_alpha(.5)
                else:
                    pass

                # facecolor ...
                if i < 4: # regular subplots always have face color
                    ax.set_facecolor(face_color)
                elif di == 'solver message': # solver message subplot always have face color
                    ax.set_facecolor(face_color)
                elif not judge1: # only have facecolor if it is not the last iteration.
                    ax.set_facecolor(face_color)
                else:
                    pass

                # ... further things.
                if judge_sd:
                    pass # May be we wanna some special sign when iteration terminated by the solver

                # ...

            # .. subplots done ...

            super_title = "mifem.MPI ITERATIONS \n> {}/{} <".format(self._computed_steps_,self._max_steps_)
            st_fontsize = 100

            alpha = (3.1 * x_len / 17.4) / y_len
            beta = 0.6 + (0.3-0.6)*(r_num_subplots-3)/2
            if judge1:
                fig.suptitle(super_title, fontsize=st_fontsize, backgroundcolor='seagreen', y=1 + alpha*beta)
            else:
                fig.suptitle(super_title, fontsize=st_fontsize, backgroundcolor='tomato', y=1 + alpha*beta)

            # noinspection PyBroadException
            try: # in case saving fail
                plt.savefig('MPI_IGR_{}.png'.format(
                    self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp), dpi=225,
                    bbox_inches='tight', facecolor='honeydew')
            except:
                pass
            plt.close(fig)
            self.___last_graph_save_time___ = time()
        else:
            pass

    def DO_send_email_to_users(self):
        """"""
        if not whether_internet_connected(): return

        # last step save and has cost long enough
        judge1 = (self._computed_steps_ == self._max_steps_ or self._iterator_.shut_down) and \
                  (self._total_cost_ > self.___email_report_time___ or self._ever_do_email_report_)
        # intermediate save: has cost a certain time or it will cost long time, so we do a first report
        judge2 = ((self._current_time_ - self.___last_email_sent_time___) > self.___email_report_time___) or \
                 self._do_first_email_warning_report_
        # to judge special iteration stop: shut down by the solver.
        judge_sd = self._iterator_.shut_down

        if self._do_first_email_warning_report_:
            self._do_first_email_warning_report_ = False
        else:
            pass

        if judge1 or judge2:

            if not self._ever_do_email_report_: self._ever_do_email_report_ = True

            if self.IS_open: return self.___PRIVATE_generate_open_email_report___()

            # noinspection PyBroadException
            try:
                html = self.summary_html
                indices = self.___PRIVATE_select_reasonable_amount_of_data___(20, last_num=5)
                RDF = self._iterator_.RDF.iloc[indices]
                RDF.to_html('{}_temp_html.html'.format(
                    self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp))
                rdf = codecs.open("{}_temp_html.html".format(
                    self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp), 'r')
                html += rdf.read()
                rdf.close()
                os.remove("{}_temp_html.html".format(
                    self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp))
                html += """
                </body>
                </html>
                """
                if judge1:
                    subject = f'mifem.MPI Completion Report {self._computed_steps_}/{self._max_steps_}'
                    header = f"""
                    <html>
                    <body>

                    <h1 style="background-color:DarkGreen;"> mifem.MPI [{self._iterator_.__class__.__name__}] 
                    HTML completion report {self._computed_steps_}/{self._max_steps_}</h1>"""
                else:
                    subject = f'mifem.MPI Processing Report {self._computed_steps_}/{self._max_steps_}'
                    header = f"""
                    <html>
                    <body>

                    <h1 style="background-color:DarkRed;"> mifem.MPI [{self._iterator_.__class__.__name__}]
                    HTML processing report {self._computed_steps_}/{self._max_steps_}</h1>"""

                if judge_sd:
                    pass #May be we wanna some special sign when iteration terminated by the solver

                html = header + html
                recent_IGR = 'MPI_IGR_{}.png'.format(
                    self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp)

                # noinspection PyBroadException
                try: # in case attaching picture fail
                    if os.path.isfile(recent_IGR):
                        SendAdminAnHTMLEmail()(html, subject=subject, attachment_images=recent_IGR)
                    else:
                        SendAdminAnHTMLEmail()(html, subject=subject)
                except:
                    SendAdminAnHTMLEmail()(html, subject=subject)

            except:  # somehow, the email sending is wrong, we then just ignore it.
                # noinspection PyBroadException
                try:
                    message = "\n<{}> is running.\n".format(self._iterator_.standard_properties.name)
                    message += "Iterations: {}/{}.\n".format(self._computed_steps_, self._max_steps_)
                    message += "Total cost: {}.\n".format(MyTimer.seconds2dhmsm(self._total_cost_))
                    message += "Each iteration costs: {}.\n".format(
                        MyTimer.seconds2dhmsm(self._average_each_run_cost_))
                    message += "Estimated remaining time: {}.\n\n".format(
                        MyTimer.seconds2dhmsm(self._estimated_remaining_time_))
                    message += "Start at :: {}.\n".format(self._str_started_time_)
                    if self._estimated_remaining_time_ <= 0.1:
                        EET = 'NOW'
                    else:
                        EET = str(self._estimated_end_time_)[:19]
                    message += "E end at:: {}.\n\n".format(EET)

                    for solver_message in self._iterator_.message:
                        message += str(solver_message) + '\n'

                    SendAdminAnEmail()(message)
                except:
                    # noinspection PyBroadException
                    try:
                        SendAdminAnEmail()("Iterations are running; sending message failed.")
                    except:
                        # all attempts failed.
                        pass # we just pass, forget about the email report.

            self.___last_email_sent_time___ = time()
        else:
            pass


    @staticmethod
    def ___DO_filter_extreme_time___(times):
        """
        We remove some extreme values to make the iteration time plot to be of higher resolution.

        Note that, after we have removed some extreme value, the iteration time plot may look
        very weird. For example, the average iteration time may be greater than all iteration times.
        This is okay, since we may have removed a huge iteration time. This should disappear after
        we have a large amount of iterations.

        :param times:
        :return:
        """
        valid_time = times[~np.isnan(times)]
        avrg = np.mean(valid_time)
        maxt = np.max(valid_time)

        if maxt > 2 * avrg and len(valid_time[valid_time>0.9*maxt]) == 1:
            # there is only one huge value. This happens in some, for example, TGV cases.
            TIME = np.nan_to_num(times)
            max_ind = np.argwhere(TIME > 0.9*maxt)[:,0]
            times[max_ind] = np.nan
        else:
            pass

        return times

