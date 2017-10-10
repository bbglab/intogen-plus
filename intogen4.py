import csv
import gzip
import os
import logging
import click
import numpy as np

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from bgparsers import selector, readers
from collections import Counter, defaultdict
from intervaltree import IntervalTree

from tasks.oncodriveclust import OncodriveClustTask
from tasks.oncodrivefml import OncodriveFmlTask
from tasks.oncodriveomega import OncodriveOmegaTask
from tasks.hotmaps import HotmapsTask
from tasks.vep import VepTask
from tasks.mutsigcv import MutsigCvTask
from tasks.schulze import SchulzeTask


TASKS = {t.KEY: t for t in [
    VepTask,
    OncodriveFmlTask,
    OncodriveOmegaTask,
    OncodriveClustTask,
    HotmapsTask,
    MutsigCvTask,
    SchulzeTask
]}


CONFIG = {
    VepTask.KEY: {'conda_env': os.environ['CONDA_ENV'], 'vep_data': os.environ['VEP_DATA']},
    OncodriveFmlTask.KEY: {'config_file': os.environ['ONCODRIVEFML_CONF']},
    OncodriveOmegaTask.KEY: {'config_file': os.environ['ONCODRIVEOMEGA_CONF']},
    HotmapsTask.KEY: {'conda_env': 'hotmaps', 'method_folder': os.environ['HOTMAPS_METHOD']}
}

logger = logging.getLogger('intogen')
LOG_FILE = 'intogen.log'
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

###############33


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


def preprocess_filtering(file, annotations=None, extra=False):

    # Minimum cutoff
    MIN_CUTOFF = 1000
    CHROMOSOMES = set(list(range(1, 23)) + ['X', 'Y'])

    # Find Hypermutators Samples
    sample_muts = Counter(
        [m['SAMPLE'] for m in readers.variants(file, annotations=annotations, extra=extra) if m['ALT_TYPE']=='snp']
    )

    vals = list(sample_muts.values())
    if len(vals) == 0:
        return

    iqr = np.subtract(*np.percentile(vals, [75, 25]))
    q3 = np.percentile(vals, 75)
    cutoff = max(MIN_CUTOFF, (q3 + 1.5 * iqr))
    hypermutators = set([k for k, v in sample_muts.items() if v > cutoff])
    if len(hypermutators) > 0:
        logger.info("[QC] {} HYPERMUTATORS at {}:  {}".format(file, cutoff, ", ".join(["{} = {}".format(h, sample_muts[h]) for h in hypermutators])))

    # Load coverage regions tree
    regions_file = os.environ['COVERAGE_REGIONS']
    coverage_tree = defaultdict(IntervalTree)
    with gzip.open(regions_file, 'rt') as fd:
        reader = csv.reader(fd, delimiter='\t')
        for i, r in enumerate(reader, start=1):
            coverage_tree[r[0]][int(r[1]):(int(r[2]) + 1)] = i

    # Read variants
    for v in readers.variants(file, annotations=annotations, extra=extra):

        # Skip hypermutators
        if v['SAMPLE'] in hypermutators:
            continue

        if v['CHROMOSOME'] not in CHROMOSOMES:
            continue

        if v['CHROMOSOME'] in coverage_tree:
            if len(coverage_tree[v['CHROMOSOME']][v['POSITION']]) == 0:
                logger.info("[QC] {} LOW COVERAGE: {} at {}:{}".format(file, v['SAMPLE'], v['CHROMOSOME'], v['POSITION']))
                continue
        yield v

######################


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option('--debug', is_flag=True, help='Enable debugging')
@click.version_option()
def cli(debug):
    if debug:
        fmt = logging.Formatter('%(asctime)s %(message)s', datefmt='%H:%M:%S')
        fh = logging.FileHandler(LOG_FILE, 'w')
        fh.setLevel(logging.DEBUG if debug else logging.INFO)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
        logger.setLevel(logging.DEBUG)
        logger.debug('Debug mode enabled')


@click.command(short_help='Create tasks input files')
@click.option('--input', '-i', required=True, help="Input file or folder", type=click.Path())
@click.option('--output', '-o', required=True, help="Output folder")
@click.option('--groupby', '-g', default="DATASET", type=str, help="Input data group by field")
@click.option('--cores', default=os.cpu_count(), type=int, help="Maximum groups to process in parallel")
@click.argument('tasks', nargs=-1)
def preprocess(input, output, groupby, cores, tasks):
    groups = selector.groupby(input, by=groupby)
    groups = list(groups)
    tasks = [TASKS[t](output, CONFIG) for t in tasks]
    prepare_tasks(groups, preprocess_filtering, tasks, cores=cores)


@click.command(short_help='Create tasks input files')
@click.option('--input', '-i', required=True, help="Input file or folder", type=click.Path())
@click.option('--output', '-o', required=True, help="Output folder")
@click.argument('tasks', nargs=-1)
def read(input, output, tasks):
    tasks = [TASKS[t](output, CONFIG) for t in tasks]
    group_key = os.path.basename(input).split('.')[0]
    prepare_tasks([(group_key, input)], read_file, tasks, cores=1)


@click.command(short_help='Run a task')
@click.option('--output', '-o', default="output", type=click.Path(), help="Output folder")
@click.argument('task', type=str)
@click.argument('key', type=str)
def run(output, task, key):
    task = TASKS[task](output, CONFIG)
    task.init(key)
    task.run()


cli.add_command(preprocess)
cli.add_command(read)
cli.add_command(run)


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    cli()
