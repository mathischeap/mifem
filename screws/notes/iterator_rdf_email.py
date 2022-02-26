# def DO_send_email_to_users(self):
#     """"""
#     if not whether_internet_connected(): return
#
#     # last step save and has cost long enough
#     judge1 = (self._computed_steps_ == self._max_steps_ or self._iterator_.shut_down) and \
#               (self._total_cost_ > self.___email_report_time___ or self._ever_do_email_report_)
#     # intermediate save: has cost a certain time or it will cost long time, so we do a first report
#     judge2 = ((self._current_time_ - self.___last_email_sent_time___) > self.___email_report_time___) or \
#              self._do_first_email_warning_report_
#     # to judge special iteration stop: shut down by the solver.
#     judge_sd = self._iterator_.shut_down
#
#     if self._do_first_email_warning_report_:
#         self._do_first_email_warning_report_ = False
#     else:
#         pass
#
#     if judge1 or judge2:
#
#         if not self._ever_do_email_report_: self._ever_do_email_report_ = True
#
#         if self.IS_open: return self.___PRIVATE_generate_open_email_report___()
#
#         # noinspection PyBroadException
#         try:
#             html = self.summary_html
#             indices = self.___PRIVATE_select_reasonable_amount_of_data___(20, last_num=5)
#             RDF = self._iterator_.RDF.iloc[indices]
#             RDF.to_html('{}_temp_html.html'.format(
#                 self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp))
#             rdf = codecs.open("{}_temp_html.html".format(
#                 self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp), 'r')
#             html += rdf.read()
#             rdf.close()
#             os.remove("{}_temp_html.html".format(
#                 self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp))
#             html += """
#             </body>
#             </html>
#             """
#             if judge1:
#                 subject = f'mifem.MPI Completion Report {self._computed_steps_}/{self._max_steps_}'
#                 header = f"""
#                 <html>
#                 <body>
#
#                 <h1 style="background-color:DarkGreen;"> mifem.MPI [{self._iterator_.__class__.__name__}]
#                 HTML completion report {self._computed_steps_}/{self._max_steps_}</h1>"""
#             else:
#                 subject = f'mifem.MPI Processing Report {self._computed_steps_}/{self._max_steps_}'
#                 header = f"""
#                 <html>
#                 <body>
#
#                 <h1 style="background-color:DarkRed;"> mifem.MPI [{self._iterator_.__class__.__name__}]
#                 HTML processing report {self._computed_steps_}/{self._max_steps_}</h1>"""
#
#             if judge_sd:
#                 #TODO: May be we wanna some special sign when iteration terminated by the solver
#                 pass
#             else:
#                 pass
#
#             html = header + html
#             recent_IGR = 'MPI_IGR_{}.png'.format(
#                 self._iterator_.standard_properties.name+'-'+self._iterator_.standard_properties.stamp)
#
#             # noinspection PyBroadException
#             try: # in case attaching picture fail
#                 if os.path.isfile(recent_IGR):
#                     SendAdminAnHTMLEmail(usErs)(html, subject=subject, attachment_images=recent_IGR)
#                 else:
#                     SendAdminAnHTMLEmail(usErs)(html, subject=subject)
#             except:
#                 SendAdminAnHTMLEmail(usErs)(html, subject=subject)
#
#         except:  # somehow, the email sending is wrong, we then just ignore it.
#             # noinspection PyBroadException
#             try:
#                 message = "\n<{}> is running.\n".format(self._iterator_.standard_properties.name)
#                 message += "Iterations: {}/{}.\n".format(self._computed_steps_, self._max_steps_)
#                 message += "Total cost: {}.\n".format(MyTimer.seconds2dhmsm(self._total_cost_))
#                 message += "Each iteration costs: {}.\n".format(
#                     MyTimer.seconds2dhmsm(self._average_each_run_cost_))
#                 message += "Estimated remaining time: {}.\n\n".format(
#                     MyTimer.seconds2dhmsm(self._estimated_remaining_time_))
#                 message += "Start at :: {}.\n".format(self._str_started_time_)
#                 if self._estimated_remaining_time_ <= 0.1:
#                     EET = 'NOW'
#                 else:
#                     EET = str(self._estimated_end_time_)[:19]
#                 message += "E end at:: {}.\n\n".format(EET)
#
#                 for solver_message in self._iterator_.message:
#                     message += str(solver_message) + '\n'
#
#                 SendAdminAnEmail(message)()
#             except:
#                 # noinspection PyBroadException
#                 try:
#                     SendAdminAnEmail("Iterations are running; sending message failed.")()
#                 except:
#                     # all attempts failed.
#                     pass # we just pass, forget about the email report.
#                     #TODO: Maybe we wanna do something when all email-sending attempts failed.
#
#         self.___last_email_sent_time___ = time()
#     else:
#         pass