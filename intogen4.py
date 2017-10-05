import os
import logging
import click

from bgparsers import selector, readers

from pipeline import read_file, prepare_tasks, run_task
from tasks.oncodriveclust import OncodriveClustTask
from tasks.oncodrivefml import OncodriveFmlTask
from tasks.oncodriveomega import OncodriveOmegaTask
from tasks.hotmaps import HotmapsTask
from tasks.vep import VepTask


TASKS = {t.KEY: t for t in [
    VepTask,
    OncodriveFmlTask,
    OncodriveOmegaTask,
    OncodriveClustTask,
    HotmapsTask
]}


CONDA_ENV = 'intogen2017'
VEP_DATA = '/home/jordeu/workspace/intogen/intogen-plus/datasets/vep'
ONCODRIVEFML_CONF = '/home/jordeu/workspace/intogen/intogen-plus/datasets/oncodrivefml/oncodrivefml.conf'
ONCODRIVEOMEGA_CONF = '/home/jordeu/workspace/intogen/intogen-plus/datasets/oncodriveomega/oncodriveomega.conf'

CONFIG = {
    VepTask.KEY: {'conda_env': CONDA_ENV, 'vep_data': VEP_DATA},
    OncodriveFmlTask.KEY: {'config_file': ONCODRIVEFML_CONF},
    OncodriveOmegaTask.KEY: {'config_file': ONCODRIVEOMEGA_CONF},
    HotmapsTask.KEY: {'conda_env': 'hotmaps', 'method_folder': '/home/jordeu/workspace/intogen/intogen-plus/methods/hotmaps'}
}

logger = logging.getLogger('intogen')
LOG_FILE = 'intogen.log'
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


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
@click.option('--groupby', '-g', default="TISSUE", type=str, help="Input data group by field")
@click.option('--cores', default=os.cpu_count(), type=int, help="Maximum groups to process in parallel")
@click.argument('tasks', nargs=-1)
def preprocess(input, output, groupby, cores, tasks):
    groups = selector.groupby(input, by=groupby)
    groups = list(groups)
    tasks = [TASKS[t](output, CONFIG) for t in tasks]
    prepare_tasks(groups, readers.variants, tasks, msg="Preprocessing", cores=cores)


@click.command(short_help='Create tasks input files')
@click.option('--input', '-i', required=True, help="Input file or folder", type=click.Path())
@click.option('--output', '-o', required=True, help="Output folder")
@click.argument('tasks', nargs=-1)
def read(input, output, tasks):
    tasks = [TASKS[t](output, CONFIG) for t in tasks]
    group_key = os.path.basename(input).split('.')[0]
    prepare_tasks([(group_key, input)], read_file, tasks, msg="Reading", cores=1)


@click.command(short_help='Run a task')
@click.option('--output', '-o', default="output", type=click.Path(), help="Output folder")
@click.argument('task', type=str)
@click.argument('key', type=str)
def run(output, task, key):
    task = TASKS[task](output, CONFIG)
    task.init(key)
    run_task(task)


cli.add_command(preprocess)
cli.add_command(read)
cli.add_command(run)


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    cli()
