
import os

from configobj import ConfigObj


CODE_DIR = os.path.dirname(os.path.abspath(__file__))
configs_file = os.path.join(CODE_DIR, 'intogen_qc.cfg')
CONF = ConfigObj(configs_file)

METHODS = list(CONF.keys())

REGIONS = os.path.join(os.environ['INTOGEN_DATASETS'], 'regions', 'cds.regions.gz')
