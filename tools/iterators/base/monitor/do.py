# -*- coding: utf-8 -*-
import numpy as np
from screws.freeze.main import FrozenOnly
import codecs
from screws.emails.plain import SendAdminAnHTMLEmail, whether_internet_connected, SendAdminAnEmail
import matplotlib, psutil, platform

import matplotlib.pyplot as plt
from matplotlib import cm
import os
from time import time, sleep
import datetime
import socket
from screws.miscellaneous.timer import MyTimer


class IteratorMonitorDo(FrozenOnly):
    """"""
    def __init__(self, monitor):
        """"""
        self._monitor_ = monitor
        self._do_first_auto_save_ = True
        memory = psutil.virtual_memory()
        self._mem_total_ = str(round(memory.total / 1024 / 1024/ 1024, 2))
        self._cpu_count_ = psutil.cpu_count(logical=False)
        self._logical_cpu_count_ = psutil.cpu_count()
        self._platform_ = platform.platform()
        self._system_ = platform.system()
        self._processor_ = platform.processor()
        self._architecture_ = platform.architecture()
        self._python_version_ = platform.python_version()
        # self._uname_ = platform.uname()

        self._freeze_self_()


    def ___PRIVATE_select_reasonable_amount_of_data___(self, max_num, last_num=1):
        """
        To report RDF, we do not report all, we make a selection.

        :param max_num: we at most selection this amount of data
        :param last_num: we must select the last ``last_num`` rows of RDF.
        """
        assert max_num >= 10
        assert 1 <= last_num < max_num
        all_data_num = len(self._monitor_._iterator_.RDF)
        if max_num < all_data_num:
            indices = list(np.linspace(0, all_data_num-last_num, max_num+1-last_num).astype(int))
            indices.extend([all_data_num+1-last_num+i for i in range(last_num-1)])
        else:
            indices = [i for i in range(all_data_num)]
        return indices

    def update(self):
        """Update monitor after every iteration."""
        monitor = self._monitor_
        monitor._computed_steps_ += 1
        if monitor._current_time_ is None:
            monitor._current_time_ = monitor._ft_firstRun_ # to make the first iteration time correct.
        monitor._last_run_cost_ = time() - monitor._current_time_
        monitor._current_time_ = time()
        monitor._total_cost_ = monitor._current_time_ - monitor._ft_start_time_

        if self._monitor_.IS.open:

            raise NotImplementedError()

        else:

            if monitor._last_run_cost_ > 0.25:  # we only consider cost long enough iteration as effective.
                monitor._effective_run_cost_ += monitor._last_run_cost_
                monitor._effective_run_num_ += 1
                monitor._average_each_run_cost_ = monitor._effective_run_cost_ / monitor._effective_run_num_
                monitor._estimated_remaining_time_ = monitor._average_each_run_cost_ * (
                        monitor._max_steps_ - monitor._computed_steps_
                )
                if monitor._estimated_remaining_time_ > 3600*24*99:
                    monitor._estimated_remaining_time_ = 3600*24*99 + 23*3600 + 59*60 + 59.999
                monitor._estimated_end_time_ = datetime.datetime.now() + datetime.timedelta(
                    seconds=monitor._estimated_remaining_time_)

                if monitor._do_first_assessment_: # will see first 5 iterations to see the estimated time.
                    monitor._do_first_assessment_counter_ += 1
                    monitor._do_first_assessment_estimate_ += monitor._estimated_remaining_time_
                    est = monitor._do_first_assessment_estimate_ / monitor._do_first_assessment_counter_
                    if est > monitor.___email_report_time___:
                        monitor._do_first_email_warning_report_ = True
                    if est > monitor.___graph_report_time___:
                        monitor._do_first_graph_warning_report_ = True
                    if (monitor._do_first_email_warning_report_ and monitor._do_first_graph_warning_report_) or \
                        monitor._do_first_assessment_counter_ >= 5:
                        monitor._do_first_assessment_ = False
                        del monitor._do_first_assessment_counter_, monitor._do_first_assessment_estimate_

                monitor._times_ = np.append(monitor._times_, monitor._last_run_cost_)
            else:
                # this is an interesting choice. No clue why I did this
                monitor._times_ = np.append(monitor._times_, np.nan)

            monitor._TIMES_ = np.append(monitor._TIMES_, monitor._total_cost_)

            if (monitor._max_steps_ == monitor._computed_steps_) or monitor._iterator_.shut_down:
                # make sure we have correct _estimated_remaining_time_ when iterator is done
                monitor._estimated_remaining_time_ = 0
                monitor._estimated_end_time_ = datetime.datetime.now()





    # noinspection PyBroadException
    def auto_save(self):
        monitor = self._monitor_

        if monitor.RDF_filename is not None:

            if monitor._last_auto_save_time_ is None:
                monitor._last_auto_save_time_ = monitor._ft_firstRun_

            gap_time = time() - monitor._last_auto_save_time_


            if monitor.auto_save_frequency is True:

                if gap_time > monitor.___auto_save_time___ or self._do_first_auto_save_:
                    try:
                        # if PermissionError, we do not stop the iteration
                        monitor._iterator_.RDF.to_csv(monitor.RDF_filename, header=True)
                    except: # wait 3 seconds
                        sleep(3)
                        try: # try once more
                            monitor._iterator_.RDF.to_csv(monitor.RDF_filename, header=True)
                        except: # just skip it
                            pass
                    monitor._last_auto_save_time_ = time()
                    self._do_first_auto_save_ = False

                else:
                    pass


            elif monitor.auto_save_frequency > 0:

                if ((monitor._computed_steps_ % monitor.auto_save_frequency == 0) and
                    gap_time > 0.1*monitor.___auto_save_time___) or \
                    gap_time > monitor.___auto_save_time___ or \
                    self._do_first_auto_save_: # we always have a look at the first iteration.
                    # important. When read from csv, it is very fast, so we do not save.
                    try: # if PermissionError, we do not stop the iteration
                        monitor._iterator_.RDF.to_csv(monitor.RDF_filename, header=True)
                    except: # wait 3 seconds
                        sleep(3)
                        try: # try once moe
                            monitor._iterator_.RDF.to_csv(monitor.RDF_filename, header=True)
                        except: # just skip it
                            pass
                    monitor._last_auto_save_time_ = time()
                    self._do_first_auto_save_ = False

                else:
                    pass


            elif self._monitor_._real_time_monitor_:

                try:
                    # if PermissionError, we do not stop the iteration
                    monitor._iterator_.RDF.to_csv(monitor.RDF_filename, header=True)
                except: # wait 3 seconds
                    sleep(3)
                    try: # try once more
                        monitor._iterator_.RDF.to_csv(monitor.RDF_filename, header=True)
                    except: # just skip it
                        pass
                monitor._last_auto_save_time_ = time()


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




    def generate_graph_report(self):
        """"""
        monitor = self._monitor_

        if monitor.IS.open: return self.___PRIVATE_generate_open_graph_report___()

        # do A LAST REPORT BEFORE STOP:
        judge1 = \
            ((monitor._computed_steps_ == monitor._max_steps_) or monitor._iterator_.shut_down) and \
            (monitor._total_cost_ > monitor.___graph_report_time___ or monitor._ever_do_graph_report_)

        # intermediate save: has cost a certain time, or it will cost long time, so we do a first report
        judge2 = ((monitor._current_time_ - monitor.___last_graph_save_time___) > monitor.___graph_report_time___) or \
            monitor._do_first_graph_warning_report_

        # real time monitoring: each single time step, we have a look
        judge3 = self._monitor_._real_time_monitor_

        # to judge special iteration stop: shut down by the solver.
        judge_sd = monitor._iterator_.shut_down

        if monitor._do_first_graph_warning_report_: monitor._do_first_graph_warning_report_ = False

        if judge1 or judge2 or judge3: # now we need to do the reporting.

            if not monitor._ever_do_graph_report_: monitor._ever_do_graph_report_ = True

            save_time = MyTimer.current_time()[1:-1]
            indices = self.___PRIVATE_select_reasonable_amount_of_data___(1000, last_num=100)
            RDF = monitor._iterator_.RDF.iloc[indices]

            matplotlib.use('Agg') # make sure we use the right backend.

            plt.rc('text', usetex=False)

            num_subplots = RDF.shape[1] + 4
            # We plot 4 extra: 't iteration', 't accumulation', solver message, and machine load

            colors = cm.get_cmap('Dark2', num_subplots-6)
            r_num_subplots = int(np.ceil(num_subplots/2))
            x_len, y_len = 18, 4.5*r_num_subplots
            fig = plt.figure(figsize=(x_len, y_len))
            plot_keys = list(RDF.columns)
            plot_keys = plot_keys[:2] + \
                        ['t iteration', 't accumulation', 'solver message', 'machine load'] + \
                        plot_keys[2:]

            # subplots ...
            for i, di in enumerate(plot_keys):
                ylabel_backgroundcolor = 'paleturquoise'
                face_color = 'aliceblue'
                ylabel = di.replace('_', '-')
                m = int(i/2)
                n = i % 2
                # noinspection PyUnboundLocalVariable
                ax = plt.subplot2grid((r_num_subplots, 2),(m, n))
                if di == 'machine load':
                    mem_percent = psutil.virtual_memory().percent
                    cpu_load = self._monitor_._iterator_.___cpu_load___
                    face_color = "whitesmoke"
                    plt.axis([0, 10, 0, 10])
                    ax.axes.get_xaxis().set_visible(False)
                    ax.axes.get_yaxis().set_visible(False)
                    plt.text(0.1, 9.5, 'MACHINE LOAD:',
                             color= 'deepskyblue', fontsize=18, style='normal',
                             ha='left', va='top', wrap=True)

                    TEXT = f"MEM: {mem_percent}% of {self._mem_total_}G used.\n" \
                           f" CPU: {cpu_load}% of {self._cpu_count_} (physical) = " \
                                f"{self._logical_cpu_count_} (logical) processors used.\n\n" \
                           f"SYSTEM: {self._system_}\n" \
                           f"{self._platform_}\n" \
                           f"{self._architecture_}\n" \
                           f"{self._processor_}\n\n" \
                           f"PYTHON version: {self._python_version_}"

                    plt.text(0.1, 8, TEXT, color= 'navy', fontsize=13,
                             ha='left', va='top', wrap=True)

                elif di == 'solver message':
                    face_color = "whitesmoke"
                    plt.axis([0, 10, 0, 10])
                    ax.axes.get_xaxis().set_visible(False)
                    ax.axes.get_yaxis().set_visible(False)
                    plt.text(0.1, 9.5, 'SOLVER MESSAGE:', color= 'deepskyblue', fontsize=18, style='normal',
                             ha='left', va='top', wrap=True)
                    message = ''
                    for M in monitor._iterator_.message:
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

                    AC = self._monitor_._iterator_.___assembling_cost___
                    if AC > 0:
                        ratio = round(AC * 100 / monitor._last_run_cost_, 2)
                        AC = round(AC, 2)
                        ITC = round(monitor._last_run_cost_, 2)
                        message += f'Assembling cost: {AC} ({ratio}% of total iteration cost: {ITC})\n\n'

                    if monitor._iterator_.shut_down:
                        message += '>>> SHUT-DOWN <<<'

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
                    t1 = f'* {monitor._iterator_.standard_properties.name}'
                    plt.text(-2, 11.3, t1, color= 'darkorchid', fontsize=26, style='normal', ha='left',
                             va='top', wrap=True)
                    sITC = MyTimer.seconds2hms(monitor._average_each_run_cost_)
                    sPPT = MyTimer.seconds2hms(monitor._preparation_time_)
                    sTTC = MyTimer.seconds2dhmsm(monitor._total_cost_) .split('.')[0]  + ']'
                    sLIC = MyTimer.seconds2hms(monitor._last_run_cost_)
                    sERT = MyTimer.seconds2dhmsm(monitor._estimated_remaining_time_).split('.')[0]  + ']'
                    t2 = 'ITC: ' + sITC  + '      LIC: '    + sLIC + '\n'
                    t2 += 'TTC: ' + sTTC + '  Preparation:' + sPPT + ' \n'
                    t2 += 'ERT: ' + sERT + '\n'
                    percentage = int(10000*(monitor._computed_steps_/monitor._max_steps_)) / 100
                    t2 += f'Iterations done: {monitor._computed_steps_}/{monitor._max_steps_} ~ {percentage}%\n'
                    t2 += f'Iterator type: {monitor._iterator_.__class__.__name__}\n'
                    plt.text(-0.5, 9.5, t2, color= 'darkblue', fontsize=22, style='normal', ha='left',
                             va='top', wrap=True)
                    t3 = 'Graph saved at:: ' + save_time + '\n'
                    t3 += 'Start at :: ' + monitor._str_started_time_ + '\n'
                    if monitor._estimated_remaining_time_ == 0:
                        t3 += 'E end at:: NOW'
                    else:
                        t3 += 'E end at:: ' + str(monitor._estimated_end_time_)[:19]
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

                    from_t = monitor._iterator_.RDF['t'][indices[0]]
                    plt.text(0, ylim[0] + (ylim[1] - ylim[0])*0.2, f"compute from: t=%.3f s."%from_t, color= 'k',
                             fontsize=22, style='normal', ha='left',
                             va='bottom', wrap=True)

                    till_t = monitor._iterator_.RDF['t'][indices[-1]]
                    plt.text(0, ylim[0] + (ylim[1] - ylim[0])*0.1, f"compute till: t=%.3f s."%till_t, color= 'k',
                             fontsize=22, style='normal', ha='left',
                             va='bottom', wrap=True)

                elif di == 't iteration':
                    ylabel_backgroundcolor = 'greenyellow'
                    face_color = 'snow'
                    itime = monitor._times_[indices]
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
                        average = monitor._average_each_run_cost_
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
                    TIME = monitor._TIMES_[indices]
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
                        sERT = MyTimer.seconds2dhmsm(monitor._estimated_remaining_time_).split('.')[0] + ']'
                        y_position = 0.9 * TIME[-1]
                        plt.text(0, y_position, 'ERT: '+ sERT, color= 'darkblue',
                                 fontsize=22, style='normal', ha='left',
                                 va='top', wrap=True)
                        # noinspection PyUnboundLocalVariable
                        if v10_ratio is not None:
                            # noinspection PyUnboundLocalVariable
                            v10_ERT_seconds = monitor._estimated_remaining_time_ * v10_ratio
                            vERT = MyTimer.seconds2dhmsm(v10_ERT_seconds).split('.')[0] + ']'
                            y_position = 0.75 * TIME[-1]
                            plt.text(0, y_position, 'V10: '+ vERT, color= 'purple',
                                     fontsize=22, style='normal', ha='left',
                                     va='top', wrap=True)

                else:
                    plt.plot(RDF['t'], RDF[di], color=colors(i-6), linewidth=1.5)

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
                elif di == 'machine load': # no ITERATING for machine load subplot
                    pass
                elif not judge1: # if not the last iteration, texture it.
                    the_text = 'ITERATING'
                    text = ax.text(0.1, 0.5, the_text, fontsize=65, color='gray',
                                   horizontalalignment='left',
                                   verticalalignment='center',
                                   transform=ax.transAxes)
                    text.set_alpha(.2)
                else:
                    pass

                # facecolor ...
                if i < 4: # regular subplots always have face color
                    ax.set_facecolor(face_color)
                elif di == 'solver message': # solver message subplot always have face color
                    ax.set_facecolor(face_color)
                elif di == 'machine load': # machine load subplot always have face color
                    ax.set_facecolor(face_color)
                elif not judge1: # only have facecolor if it is not the last iteration.
                    ax.set_facecolor(face_color)
                else:
                    pass

                # ... further things.
                if judge_sd:
                    pass # Maybe we wanna some special sign when iteration terminated by the solver

                # ...

            # .. subplots done ...

            super_title = "mifem.MPI ITERATIONS \n> {}/{} <".format(monitor._computed_steps_,monitor._max_steps_)
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
                    monitor._iterator_.standard_properties.name), dpi=225,
                    bbox_inches='tight', facecolor='honeydew')
            except:
                pass
            plt.close(fig)
            monitor.___last_graph_save_time___ = time()
        else:
            pass

    def send_email_to_users(self):
        """"""
        if not whether_internet_connected(): return

        hostname = socket.gethostname()
        if hostname not in ('DT-YI-HT20', 'DESKTOP-SYSU-YiZhang'): return



        monitor = self._monitor_

        if monitor.IS.open: return self.___PRIVATE_generate_open_email_report___()




        # last step save and has cost long enough
        judge1 = (monitor._computed_steps_ == monitor._max_steps_ or monitor._iterator_.shut_down) and \
                  (monitor._total_cost_ > monitor.___email_report_time___ or monitor._ever_do_email_report_)
        # intermediate save: has cost a certain time, or it will cost long time, so we do a first report
        judge2 = ((monitor._current_time_ - monitor.___last_email_sent_time___) > monitor.___email_report_time___) or \
                 monitor._do_first_email_warning_report_
        # to judge special iteration stop: shut down by the solver.
        judge_sd = monitor._iterator_.shut_down

        if monitor._do_first_email_warning_report_:
            monitor._do_first_email_warning_report_ = False
        else:
            pass

        if judge1 or judge2:

            if not monitor._ever_do_email_report_: monitor._ever_do_email_report_ = True


            # noinspection PyBroadException
            try:
                html = monitor.summary_html
                indices = self.___PRIVATE_select_reasonable_amount_of_data___(20, last_num=5)
                RDF = monitor._iterator_.RDF.iloc[indices]
                RDF.to_html('{}_temp_html.html'.format(
                    monitor._iterator_.standard_properties.name))
                rdf = codecs.open("{}_temp_html.html".format(
                    monitor._iterator_.standard_properties.name), 'r')
                html += rdf.read()
                rdf.close()
                os.remove("{}_temp_html.html".format(
                    monitor._iterator_.standard_properties.name))
                html += """
                </body>
                </html>
                """
                if judge1:
                    subject = f'mifem.MPI Completion Report {monitor._computed_steps_}/{monitor._max_steps_}'
                    header = f"""
                    <html>
                    <body>

                    <h1 style="background-color:DarkGreen;"> mifem.MPI [{monitor._iterator_.__class__.__name__}] 
                    HTML completion report {monitor._computed_steps_}/{monitor._max_steps_}</h1>"""
                else:
                    subject = f'mifem.MPI Processing Report {monitor._computed_steps_}/{monitor._max_steps_}'
                    header = f"""
                    <html>
                    <body>

                    <h1 style="background-color:DarkRed;"> mifem.MPI [{monitor._iterator_.__class__.__name__}]
                    HTML processing report {monitor._computed_steps_}/{monitor._max_steps_}</h1>"""

                if judge_sd:
                    pass #May be we wanna some special sign when iteration terminated by the solver

                html = header + html
                recent_IGR = 'MPI_IGR_{}.png'.format(
                    monitor._iterator_.standard_properties.name)

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
                    message = "\n<{}> is running.\n".format(monitor._iterator_.standard_properties.name)
                    message += "Iterations: {}/{}.\n".format(monitor._computed_steps_, monitor._max_steps_)
                    message += "Total cost: {}.\n".format(MyTimer.seconds2dhmsm(monitor._total_cost_))
                    message += "Each iteration costs: {}.\n".format(
                        MyTimer.seconds2dhmsm(monitor._average_each_run_cost_))
                    message += "Estimated remaining time: {}.\n\n".format(
                        MyTimer.seconds2dhmsm(monitor._estimated_remaining_time_))
                    message += "Start at :: {}.\n".format(monitor._str_started_time_)
                    if monitor._estimated_remaining_time_ <= 0.1:
                        EET = 'NOW'
                    else:
                        EET = str(monitor._estimated_end_time_)[:19]
                    message += "E end at:: {}.\n\n".format(EET)

                    for solver_message in monitor._iterator_.message:
                        message += str(solver_message) + '\n'

                    SendAdminAnEmail()(message)
                except:
                    # noinspection PyBroadException
                    try:
                        SendAdminAnEmail()("Iterations are running; sending message failed.")
                    except:
                        # all attempts failed.
                        pass # we just pass, forget about the email report.

            monitor.___last_email_sent_time___ = time()
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
