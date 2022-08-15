# -*- coding: utf-8 -*-

import json
from urllib.request import urlopen


def get_my_IP_data():
    """

    :return:
    """
    url = 'http://ipinfo.io/json'
    response = urlopen(url)
    data = json.load(response)
    return data