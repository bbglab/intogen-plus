"""
This module contains code related with the configuration file (see :ref:`project configuration`).

Additionally, it includes other file realted code, specially from :mod:`bgconfig`.
"""

import logging
import os
import sys
from bgconfig import BGConfig, _file_name, _file_exists_or_die

from simReg import __logger_name__
logger = logging.getLogger(__logger_name__)

file_exists_or_die = _file_exists_or_die
file_name = _file_name


def load_configuration(config_file, override=None):
    """
    Load the configuration file and checks the format.

    Args:
        config_file: configuration file path

    Returns:
        :class:`bgconfig.BGConfig`: configuration as a :obj:`dict`

    """
    config_template = os.path.join(os.path.dirname(__file__), "oncodrivefml_v2.conf.template")

    try:
        return BGConfig(config_template, config_file=config_file, use_env_vars=False, override_values=override, unrepr=True)
    except ValueError as e:
        logger.error(e)
        sys.exit(-1)


possible_extensions = ['.gz', '.xz', '.bz2', '.tsv', '.txt']
"""
Some expected extensions
"""
special_characters = ['.', '_']
"""
Some special characters
"""
def remove_extension_and_replace_special_characters(file_path):
    """
    Modifies the name of a file by removing any extension in :attr:`possible_extensions`
    and replacing any character in :attr:`special_characters` for ``-``.

    Args:
        file_path: path to a file

    Returns:
        str: file name modified

    """
    file_name=os.path.basename(file_path)
    for ext in possible_extensions:
        file_name = file_name.replace(ext, '')
    for char in special_characters:
        file_name = file_name.replace(char, '-')
    return file_name
