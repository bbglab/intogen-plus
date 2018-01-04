"""
This module contains code related to the configuration file (see :ref:`project configuration`).
"""

import logging
import os
import sys
from bgconfig import BGConfig


def load_configuration(config_file=None):
    """
    Load the configuration file and checks the format.
    Returns:
        :class:`bgconfig.BGConfig`: configuration as a :obj:`dict`
    """
    config_template = os.path.join(os.path.dirname(__file__), "mutrate.conf.template")
    config_spec = os.path.join(os.path.dirname(__file__), "mutrate.conf.template.spec")
    try:
        config = BGConfig(config_template, config_file=config_file, config_spec=config_spec)
        return config
    except ValueError as e:
        logging.error(e)
        sys.exit(-1)