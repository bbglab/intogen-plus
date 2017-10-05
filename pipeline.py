import csv
import gzip
import logging
import signal
import sys

from os import path, kill, getpgid, killpg, makedirs
from concurrent.futures import ProcessPoolExecutor
from copy import copy
from datetime import datetime
from functools import partial
from time import sleep

from ago import human
from tqdm import tqdm
from io import TextIOWrapper

LOG = logging.getLogger('intogen')

G_STDERR = sys.stderr
G_STDOUT = sys.stdout
G_HANDLERS = logging.getLogger().handlers
G_DEFAULTS = tqdm.__init__.__defaults__


def read_file(tsv_file):
    with gzip.open(tsv_file, "rt") as fd:
        reader = csv.DictReader(fd, delimiter='\t')
        for row in reader:
            yield row


def _restore_output():
    # Restore outputs
    sys.stderr = G_STDERR
    sys.stdout = G_STDOUT
    logging.getLogger().handlers = G_HANDLERS
    tqdm.__init__.__defaults__ = G_DEFAULTS


def run_task(task):
    _restore_output()
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    LOG.info("{} [start]".format(task))
    t_ini = datetime.now()

    # Redirect any output to the log file
    log_folder = path.join(task.output_folder, "logs")
    makedirs(log_folder, exist_ok=True)
    log_file = open(path.join(log_folder, "{}.log".format(task.name)), "wt")
    sys.stderr = log_file
    sys.stdout = log_file
    logging.getLogger().handlers = [logging.StreamHandler(log_file)]
    tqdm.__init__.__defaults__ = tuple(log_file if isinstance(v, TextIOWrapper) else v for v in tqdm.__init__.__defaults__)

    # Execute the task
    error = None
    try:
        task.run()
    except BaseException as e:
        error = e

    # Report errors
    _restore_output()
    if error is not None:
        LOG.error("{} [error] - '{}'".format(task, error))
        LOG.exception(error)

    log_file.close()

    return task, t_ini - datetime.now()


def run_tasks(tasks, msg=None):

    if msg:
        LOG.info(msg)

    tasks_done = set()

    with ProcessPoolExecutor(max_workers=4) as executor:
        try:
            for task, dt in executor.map(run_task, tasks):
                LOG.info("{} [done] ({})".format(task, human(dt)))
                tasks_done.add(task)
        except:

            LOG.error("Canceling processes. Please wait.")

            processes = executor._processes
            executor.shutdown(wait=False)
            for pid, process in processes.items():
                process.terminate()

            tasks_not_done = {t for t in tasks if t not in tasks_done}
            for t in tasks_not_done:
                LOG.error("Cleaning unfinished {}".format(t))
                t.clean()

            raise


def prepare_task(reader, tasks, item):
    position, (group_key, group_data) = item

    # Initialize task name
    [t.init(group_key) for t in tasks]

    # Input start
    [t.input_start() for t in tasks]

    # Input write
    for i, mut in tqdm(enumerate(reader(group_data)), desc="Reading '{}'".rjust(20).format(group_key)):
        id = "I{:010d}".format(i)
        for t in tasks:
            t.input_write(id, mut)

    # Input close
    [t.input_end() for t in tasks]

    return [(t.KEY, group_key) for t in tasks]


def prepare_tasks(groups, reader, tasks, msg=None, cores=None):

    if msg:
        LOG.info(msg)

    func = partial(prepare_task, reader, tasks)

    all_tasks = []
    if cores > 1:
        with ProcessPoolExecutor(cores) as executor:
            for tasks in executor.map(func, enumerate(groups)):
                all_tasks += tasks
    else:
        for tasks in map(func, enumerate(groups)):
            all_tasks += tasks

    return all_tasks


def execute(groups, reader, tasks, msg=None):
    keys, tasks = prepare_tasks(groups, tasks, reader, msg=msg)
    run_tasks(tasks, msg=None)
    return keys, tasks

