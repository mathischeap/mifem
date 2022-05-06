# -*- coding: utf-8 -*-
import socket
from root.config.main import rAnk, mAster_rank


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