
import logging
from datetime import datetime

from ago import human


logger = logging.getLogger(__name__)


def executor_run(executor):
    """
    Method to call the run method

    Args:
        executor (:class:`~smregions.executor.ElementExecutor`):

    Returns:
        :meth:`~smregions.executor.ElementExecutor.run`

    """
    return executor.run()


def loop_logging(iterable, size=None, step=1):
    """
    Loop through an iterable object displaying messages
    using :func:`~logging.info` on the appropiate logger

    Args:
        iterable:
        size (int): Defaults to None.
        step (int): Defaults to 1.

    Yields:
        The iterable element

    """

    if size is None:
        size = len(iterable)

    i = 0
    start_time = datetime.now()
    for i, value in enumerate(iterable):
        if i % step == 0:
            logger.info("[%d of %d]", i+1, size)
        yield value
    logger.info("[%d of %d]", i+1, size)
    logger.debug("Time: %s", human(start_time))
