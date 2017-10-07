import csv
import gzip

from concurrent.futures import ProcessPoolExecutor
from functools import partial


def read_file(tsv_file):
    with gzip.open(tsv_file, "rt") as fd:
        reader = csv.DictReader(fd, delimiter='\t')
        for row in reader:
            yield row


def prepare_task(reader, tasks, item):
    position, (group_key, group_data) = item

    # Initialize task name
    [t.init(group_key) for t in tasks]

    # Input start
    [t.input_start() for t in tasks]

    # Input write
    for i, mut in enumerate(reader(group_data)):
        id = "I{:010d}".format(i)
        for t in tasks:
            t.input_write(id, mut)

    # Input close
    [t.input_end() for t in tasks]

    return [(t.KEY, group_key) for t in tasks]


def prepare_tasks(groups, reader, tasks, cores=None):
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



