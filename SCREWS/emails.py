


import sys
if './' not in sys.path: sys.path.append('./')
import os
from SCREWS.frozen import FrozenOnly
from SCREWS.miscellaneous import MyTimer
import socket
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from root.config import rAnk, mAster_rank


class EmailSendingError(Exception):
    """Raise when we try to access standard property ``statistic`` but ``___statistic___`` is not defined."""


class SendAdminAnEmail(FrozenOnly):
    """To send an email to the library admin."""
    def __init__(self):
        assert rAnk == mAster_rank, "Should only call it in master core."
        self._MESSAGE_ = None
        self._freeze_self_()

    def ___generate_MESSAGE_2b_sent___(self, message, key_code):
        local_IP = socket.gethostbyname(socket.gethostname())
        local_machine_name = socket.gethostname()
        MESSAGE = 'Dear admin,\n\n'
        MESSAGE += 'A message from <' + local_machine_name + '@' + local_IP + '>:\n'
        MESSAGE += 'KEY-CODE: ' + key_code + '\n\n'
        MESSAGE += "“" + message + "”\n\n"
        MESSAGE += 'Best regards,\n'
        MESSAGE += 'mifem\n'
        MESSAGE += 'Hosted at: https://gitlab.com/zhangyi_aero/mifem\n'
        MESSAGE += 'Contact: zhangyi_aero@hotmail.com or y.zhang-14@tudelft.nl'
        self._MESSAGE_ = MESSAGE

    def __call__(self, message, subject=None):
        """By calling the object, we send the message."""
        assert rAnk == mAster_rank, "Should only call it in master core."
        hostname = socket.gethostname()

        # noinspection PyBroadException
        try:
            with open('root/___private_developer_code___.txt', 'r') as f:
                PDCs = f.readlines()
            key_code = PDCs[0]
            DC_abc, DC_bbb, DC_from, DC_name, DC_email1, DC_email2, DC_HS = PDCs[1:8]
        except:
            return 0

        if not whether_internet_connected():
           return 0

        in_the_wall = whether_in_the_great_fire_wall()

        # ---- connect to the serve -----------------------------------
        if in_the_wall and hostname == 'DT-YI-HT20': # edit or add condition here when I am working from new machine (only for the library holder)
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
        #===================================================================================

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




class SendAdminAnHTMLEmail(FrozenOnly):
    """To send an email to the library admin."""
    def __init__(self):
        assert rAnk == mAster_rank, "Should only call it in master core."
        hostname = socket.gethostname()

        # noinspection PyBroadException
        try:
            with open('root/___private_developer_code___.txt', 'r') as _:
                pass
        except:
            self.s = None
            self._freeze_self_()
            return

        if not whether_internet_connected():
            self.s = None

        else:

            in_the_wall = whether_in_the_great_fire_wall()
            #---- connect to the serve -----------------------------------
            if in_the_wall and hostname == 'DT-YI-HT20': # edit or add condition here when I am working from new machine (only for the library holder)
                # noinspection PyBroadException
                try:
                    with open('root/___private_developer_code___.txt', 'r') as f:
                        PDCs = f.readlines()
                    self._key_code_ = PDCs[0]
                    DC_abc, DC_bbb, DC_from, DC_name, DC_email1, DC_email2, DC_HS = PDCs[1:8]
                    # build connection: I am in the wall, so I use QQ mail. So far I use my own QQ mail.
                    host_server = DC_HS[:-1]
                    abc = DC_abc[:-1]
                    smtp = smtplib.SMTP_SSL(host_server)
                    # smtp.set_debuglevel(1) # un-comment it to show debug information.
                    smtp.ehlo(host_server)
                    smtp.login(DC_bbb[:-1], abc)
                    self.s = smtp
                    self._from_ = DC_from[:-1]

                except:
                    self.s = None
            else:
                self.s = None

            # =====================================================================

        self._names_ = list()
        self._emails_ = list()

        if self.s is not None:
            # noinspection PyUnboundLocalVariable
            self._names_ = [DC_name[:-1], DC_name[:-1]]
            # noinspection PyUnboundLocalVariable
            self._emails_ = [DC_email1[:-1], DC_email2[:-1]]
        else:
            pass

        self._freeze_self_()

    def __call__(self, HTML_message, subject=None, attachment_images=None):
        """By calling the object, we send the message."""
        assert rAnk == mAster_rank, "Should only call it in master core."

        if self.s is None:
            return 0

        # noinspection PyBroadException
        try:
            local_IP = socket.gethostbyname(socket.gethostname())
            local_machine_name = socket.gethostname()

            if attachment_images is None:
                attachment_images = list()
            if isinstance(attachment_images, str):
                attachment_images = [attachment_images, ]

            for name, email in zip(self._names_, self._emails_):
                msg = MIMEMultipart()  # create a message class
                msg['From'] = self._from_
                msg['To'] = email
                if subject is None:
                    msg['Subject'] = "mifem message (please do not reply)"
                else:
                    msg['Subject'] = subject
                MST = HTML_message + f"""
    
                <html>
                  <head></head>
                  <body>
                    <p>Best regards,<br>
                       <em> <font color="cyan">mifem</font> </em> <br>
                       Hosted at: <font color="blue">https://gitlab.com/zhangyi_aero/mifem</font> <br>
                       Contact: <font color="blue">zhangyi_aero@hotmail.com or y.zhang-14@tudelft.nl</font> <br>
                       Email composed at <font color="Crimson">{MyTimer.current_time()[1:-1]}</font> from host: < {local_machine_name} @ {local_IP} > <br>
                       KEY-CODE: {self._key_code_}
                    </p>
                  </body>
                </html>
                """
                part1 = MIMEText(MST, 'html')
                msg.attach(part1)
                for image_filename in attachment_images:
                    f = open(image_filename, 'rb')
                    img_data = f.read()
                    image = MIMEImage(img_data, name=os.path.basename(image_filename))
                    msg.attach(image)
                    f.close()
                self.s.send_message(msg)
            self.s.quit()
            return 1

        except:
            return 0









def whether_internet_connected(domain_names=None):
    """
    :param domain_names: (`default`:``None``) Domains that we use to check the internet connection.
    :type: list, tuple, None
    :return: ``True`` if we have internet connection.
    :rtype: bool
    """
    assert rAnk == mAster_rank, "Should only call it in master core."
    if domain_names is None:  # we use following default domain names.
        domain_names = ["www.qq.com", "www.baidu.com", "www.google.com", "www.microsoft.com", "www.apple.com"]
    # _____ check domain_names ...
    if isinstance(domain_names, str):
        domain_names = [domain_names, ]
    elif isinstance(domain_names, (tuple, list)):
        for dn in domain_names:
            assert isinstance(dn, str), "A domain name must be str."
    else:
        raise Exception("domain_names={} wrong.".format(domain_names))
    # ____ check internet connection ...
    connected = [False for _ in range(len(domain_names))]
    for i, dn in enumerate(domain_names):
        try:
            # see if we can resolve the host name -- tells us if there is a DNS listening
            host = socket.gethostbyname(dn)
            # connect to the host -- tells us if the host is actually reachable
            s = socket.create_connection((host, 80), 2)
            s.close()
            connected[i] = True
        except socket.error:
            connected[i] = False
        if connected[i]:
            break
    # ...
    return any(connected)


def whether_in_the_great_fire_wall():
    """Return True if we are in the great cyber fire wall."""
    assert rAnk == mAster_rank, "Should only call it in master core."
    # we use these typical walled websites to check.
    domain_names = ["www.google.com", "www.facebook.com", "www.twitter.com"]
    IN_THE_WALL = [True for _ in range(len(domain_names))]
    for i, dn in enumerate(domain_names):
        try:
            # see if we can resolve the host name -- tells us if there is a DNS listening
            host = socket.gethostbyname(dn)
            # connect to the host -- tells us if the host is actually reachable
            s = socket.create_connection((host, 80), 2)
            s.close()
            IN_THE_WALL[i] = False
        except socket.error:
            IN_THE_WALL[i] = True
    return all(IN_THE_WALL)




if __name__ == '__main__':
    # mpiexec -n 8 python SCREWS\emails.py
    if rAnk == mAster_rank:
        message = '123test'
        code = SendAdminAnEmail()(message)
        print(code)
        message = '123test'
        code = SendAdminAnHTMLEmail()(message)
        print(code)