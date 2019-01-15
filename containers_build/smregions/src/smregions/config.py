"""
This module contains code related with the configuration file.

Additionally, it includes other file related code, specially from :mod:`bgconfig`.
"""

import logging
import os
import sys

from bgconfig import BGConfig, _file_name, _file_exists_or_die


file_exists_or_die = _file_exists_or_die
file_name = _file_name


def load(config_file, override=None):
    """
    Load the configuration file and checks the format.

    Args:
        config_file: configuration file path

    Returns:
        :class:`bgconfig.BGConfig`: configuration as a :obj:`dict`

    """
    config_template = os.path.join(os.path.dirname(__file__), "smregions.conf.template")

    try:
        return BGConfig(config_template, config_file=config_file, use_env_vars=True, override_values=override, unrepr=False)
    except ValueError as e:
        logging.getLogger(__name__).error(e)
        sys.exit(-1)
