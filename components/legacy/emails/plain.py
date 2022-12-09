
# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('/')
from components.freeze.main import FrozenOnly
import socket, os
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from root.config.main import RANK, MASTER_RANK
from components.legacy.emails.html import SendAdminAnHTMLEmail
from components.legacy.internet.connected_or_not import whether_internet_connected
from components.legacy.internet.in_the_wall_or_not import whether_in_the_great_fire_wall



class SendAdminAnEmail(FrozenOnly):
    """Send emails to the library admin."""
    def __init__(self):
        assert RANK == MASTER_RANK, "Should only call it in master core."
        self._MESSAGE_ = None
        self._freeze_self_()

    def ___generate_MESSAGE_2b_sent___(self, message, key_code):
        local_IP = socket.gethostbyname(socket.gethostname())
        local_machine_name = socket.gethostname()
        MESSAGE = 'Dear admin,\n\n'
        MESSAGE += 'A message from <' + local_machine_name + '@' + local_IP + '>\n'
        MESSAGE += 'with KEY-CODE: ' + key_code[:-1] + ':' + '\n\n'
        MESSAGE += "“" + message + "”\n\n"
        MESSAGE += 'Best regards,\n'
        MESSAGE += 'mifem\n'
        MESSAGE += 'Hosted at: https://github.com/mathischeap/mifem\n'
        MESSAGE += 'Contact: zhangyi_aero@hotmail.com or zhangyi55@mail.sysu.edu.cn'
        self._MESSAGE_ = MESSAGE

    def __call__(self, message, subject=None):
        """By calling the object, we send the message."""
        assert RANK == MASTER_RANK, "Should only call it in master core."
        hostname = socket.gethostname()

        if hostname not in ('DT-YI-HT20', 'DESKTOP-SYSU-YiZhang'):
            return 0

        if not whether_internet_connected():
           return 0

        # noinspection PyBroadException
        absolute_path = os.path.dirname(__file__)

        try:
            with open(absolute_path + '/___private_developer_code___.txt', 'r') as f:
                PDCs = f.readlines()
            key_code = PDCs[0]
            DC_abc, DC_bbb, DC_from, DC_name, DC_email1, DC_email2, DC_HS = PDCs[1:8]
        except FileNotFoundError:
            return 0

        in_the_wall = whether_in_the_great_fire_wall()

        if in_the_wall:
            # noinspection PyBroadException
            try:
                # build connection: I am in the wall!
                host_server = DC_HS[:-1]
                abc = DC_abc[:-1]
                smtp = smtplib.SMTP_SSL(host_server)
                # smtp.set_debuglevel(1) # un-comment it to show debug information.
                smtp.ehlo(host_server)
                smtp.login(DC_bbb[:-1], abc)
                _from_ = DC_from[:-1]
            except:
                return 0
        else:
            return 0

        # noinspection PyBroadException
        try:
            names = (DC_name[:-1], DC_name[:-1])
            emails = (DC_email1[:-1], DC_email2[:-1])
            self.___generate_MESSAGE_2b_sent___(message, key_code)
            assert self._MESSAGE_ is not None, "I have no MESSAGE to be sent at all."

            for name, email in zip(names, emails):
                msg = MIMEMultipart()  # create a message class
                msg['From'] = _from_
                msg['To'] = email
                if subject is None:
                    msg['Subject'] = "mifem message (please do not reply)"
                else:
                    msg['Subject'] = subject
                msg.attach(MIMEText(self._MESSAGE_))
                smtp.send_message(msg)
            smtp.quit()
            return 1

        except:
            return 0









if __name__ == '__main__':
    # mpiexec -n 8 python SCREWS\emails\plain.py
    if RANK == MASTER_RANK:
        message = '123test'
        code = SendAdminAnEmail()(message)
        print(code)


        SENDER = SendAdminAnHTMLEmail
        # message = '123test'
        # code = SendAdminAnHTMLEmail()(message)
        # print(code)