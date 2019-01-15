import os
import logging
import click

from tasks.base import prepare_task, prepare_tasks
from tasks.cbase import CBaseTask
from tasks.deconstructsig import DeconstructSigTask
from tasks.dndscv import DndsCvTask
from tasks.oncodriveclustl import OncodriveClustlTask
from tasks.oncodrivefml import OncodriveFmlTask
from tasks.hotmaps import HotmapsTask
from tasks.vep import VepTask
from tasks.combination import CombinationTask
from tasks.mutrate import MutrateTask
from tasks.smregions import SmregionsTask
from filters.base import VariantsReader, TSVReader
from filters.variants import VariantsFilter
from filters.vep import VepFilter
from bgparsers import selector



TASKS = {t.KEY: t for t in [
    VepTask,
    OncodriveFmlTask,
    OncodriveClustlTask,
    HotmapsTask,
    DeconstructSigTask,
    CombinationTask,
    DndsCvTask,
    CBaseTask,
    MutrateTask,
    SmregionsTask
]}


logger = logging.getLogger('intogen')
LOG_FILE = 'intogen.log'


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--debug', is_flag=True, help='Enable debugging')
@click.version_option()
def cli(debug=False):
    if debug:
        fmt = logging.Formatter('%(asctime)s %(message)s', datefmt='%H:%M:%S')
        fh = logging.FileHandler(LOG_FILE, 'w')
        fh.setLevel(logging.DEBUG if debug else logging.INFO)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
        logger.setLevel(logging.DEBUG)
        logger.debug('Debug mode enabled')


@click.command(short_help='Create tasks input files from VARIANTS datasets')
@click.option('--input', '-i', required=True, help="Input file or folder", type=click.Path())
@click.option('--output', '-o', required=True, help="Output folder")
@click.option('--groupby', '-g', default="DATASET", type=str, help="Input data group by field")
@click.option('--cores', default=os.cpu_count(), type=int, help="Maximum groups to process in parallel")
@click.argument('tasks', nargs=-1)
def readvariants(input, output, groupby, cores, tasks):
    groups = selector.groupby(input, by=groupby)
    groups = list(groups)
    tasks = [TASKS[t](output) for t in tasks]

    reader = VariantsReader()
    for f in [VariantsFilter]:
        reader = f(reader)

    prepare_tasks(output, groups, reader, tasks, cores=cores)


@click.command(short_help='Create tasks input files from VEP output')
@click.option('--input', '-i', required=True, help="Input file or folder", type=click.Path())
@click.option('--output', '-o', required=True, help="Output folder")
@click.argument('tasks', nargs=-1)
def readvep(input, output, tasks):
    tasks = [TASKS[t](output) for t in tasks]
    group_key = os.path.basename(input).split('.')[0]
    reader = TSVReader()
    for f in [VepFilter]:
        reader = f(reader)

    prepare_tasks(output, [(group_key, input)], reader, tasks, cores=1)


@click.command(short_help='Run a task')
@click.option('--output', '-o', default="output", type=click.Path(), help="Output folder")
@click.argument('task', type=str)
@click.argument('key', type=str)
def run(output, task, key):
    task = TASKS[task](output)
    task.init(key)
    task.run()


cli.add_command(readvariants)
cli.add_command(readvep)
cli.add_command(run)


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    cli()
