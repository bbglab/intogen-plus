import logging

from os import path
from bgparsers import selector, readers

from pipeline import read_file, execute
from tasks.oncodriveclust import OncodriveClustTask
from tasks.oncodrivefml import OncodriveFmlTask
from tasks.oncodriveomega import OncodriveOmegaTask
from tasks.vep import VepTask

CONDA_ENV = '/home/jordeu/workspace/intogen/intogen-plus/envs/intogen4'
VEP_DATA = '/home/jordeu/workspace/intogen/intogen-plus/datasets/vep'
ONCODRIVEFML_CONF = '/home/jordeu/workspace/intogen/intogen-plus/datasets/oncodrivefml/oncodrivefml.conf'
ONCODRIVEOMEGA_CONF = '/home/jordeu/workspace/intogen/intogen-plus/datasets/oncodriveomega/oncodriveomega.conf'
INPUT = '/home/jordeu/workspace/intogen/intogen-plus/input/pancanatlas'
OUTPUT = '/home/jordeu/workspace/intogen/intogen-plus/output'


LOG = logging.getLogger('intogen')


def run():

    where = None
    where = {'TISSUE': {'UVM'}}

    groups = selector.groupby(INPUT, by='TISSUE', where=where)
    keys, tasks = execute(
        groups,
        readers.variants,
        [
            VepTask(output_folder=OUTPUT, conda_env=CONDA_ENV, vep_data=VEP_DATA),
            OncodriveFmlTask(output_folder=OUTPUT, config_file=ONCODRIVEFML_CONF)
        ],
        "STEP1"
    )

    keys = ['UVM']

    keys, tasks = execute(
        [(k, path.join(OUTPUT, "vep", "{}.out.gz".format(k))) for k in keys],
        read_file,
        [
            OncodriveClustTask(output_folder=OUTPUT),
            OncodriveOmegaTask(output_folder=OUTPUT, config_file=ONCODRIVEOMEGA_CONF)
        ],
        "STEP2"
    )


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    run()
