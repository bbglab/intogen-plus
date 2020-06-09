import contextlib
import gzip
import sys


@contextlib.contextmanager
def out_open(file=None, mode='w'):
    """
    Open file in write mode or send to STDOUT
    if file is None.
    Args:
        file: file path
        mode:
    Returns:
        File descriptor.
    """
    if file is None:
        fh = sys.stdout
    elif file.endswith('.gz'):
        fh = gzip.open(file, mode)
    else:
        fh = open(file, mode)
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()
