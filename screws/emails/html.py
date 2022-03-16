


from screws.freeze.main import FrozenOnly
import os
import socket
import smtplib
from email.mime.image import MIMEImage
from email.mime.text import MIMEText
from screws.miscellaneous.timer import MyTimer
from root.config.main import rAnk, mAster_rank
from email.mime.multipart import MIMEMultipart
from screws.internet.in_the_wall_or_not import whether_in_the_great_fire_wall
from screws.internet.connected_or_not import whether_internet_connected




class SendAdminAnHTMLEmail(FrozenOnly):
    """Send emails to the library admin."""

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
            # edit or add condition here when I am working from new machine (only for the library holder)
            if in_the_wall and hostname in ('DT-YI-HT20', 'DESKTOP-SYSU-YiZhang'):
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
                       Hosted at: <font color="blue">https://github.com/mathischeap/mifem</font> <br>
                       Contact: <font color="blue">zhangyi_aero@hotmail.com or zhangyi55@mail.sysu.edu.cn</font> <br>
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

