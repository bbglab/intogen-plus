import logging

from ago import human
from datetime import datetime
from tqdm import tqdm
from bgparsers import selector, readers

from tasks.oncodrivefml import OncodriveFmlTask
from tasks.vep import VepTask

CONDA_ENV = '/home/jordeu/workspace/intogen/intogen-plus/envs/intogen4'
VEP_DATA = '/home/jordeu/workspace/intogen/intogen-plus/datasets/vep'
ONCODRIVEFML_CONF = '/home/jordeu/workspace/intogen/intogen-plus/datasets/oncodrivefml/oncodrivefml.conf'
INPUT = '/home/jordeu/workspace/intogen/intogen-plus/input/pancanatlas'
OUTPUT = '/home/jordeu/workspace/intogen/intogen-plus/output'


LOG = logging.getLogger('intogen')


def step1_prepare(group):
    group_key, group_selection = group

    # Task creation
    tasks = [
        VepTask(name=group_key, output_folder=OUTPUT, conda_env=CONDA_ENV, vep_data=VEP_DATA),
        OncodriveFmlTask(name=group_key, output_folder=OUTPUT, config_file=ONCODRIVEFML_CONF)
    ]

    # Input start
    [t.input_start() for t in tasks]

    # Input write
    for i, mut in tqdm(enumerate(readers.variants(group_selection)), desc="Reading '{}'".rjust(20).format(group_key)):
        id = "V{:010d}".format(i)
        for t in tasks:
            t.input_write(id, mut)

    # Input close
    [t.input_end() for t in tasks]

    return tasks


def step2_run(task):
    LOG.info("{} [start]".format(task))
    t_ini = datetime.now()
    task.run()
    return task, datetime.now() - t_ini


def run():

    where = {'TISSUE': {'UVM', 'PCPG', 'KICH', 'TGCT'}}
    # where = {'TISSUE': {'PCPG'}}

    # STEP1. Prepare input files
    LOG.info("STEP1. Preprocessing input files")
    all_tasks = []
    for tasks in map(step1_prepare, selector.groupby(INPUT, by='TISSUE', where=where)):
        all_tasks += tasks

    # STEP2. Run tasks
    LOG.info("STEP2. Running VEP and OncodriveFML")
    for task, dt in map(step2_run, all_tasks):
        LOG.info("{} [done] ({})".format(task, human(dt)))


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    run()

