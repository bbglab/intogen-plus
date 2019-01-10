import logging
import numpy as np

import bgdata
from bgpack import packer, reader as bgpack_reader

from simReg import __logger_name__
from simReg.signature import get_ref_triplet
from simReg.walker import partitions_list

logger = logging.getLogger(__logger_name__)


class ElementExecutor:
    """
    Base class for executors that do the analysis per genomic element.
    The mutations are simulated to compute a p_value
    by comparing the functional impact of the observed mutations
    and the expected.
    The simulation parameters are taken from the configuration file.

    Args:
        element_id (str): element ID
        muts (list): list of mutations belonging to that element (see :ref:`mutations <mutations dict>` inner items)
        segments (list): list of segments belonging to that element (see :ref:`elements <elements dict>` inner items)
        signature (dict): probabilities of each mutation (see :ref:`signatures <signature dict>`)
        config (dict): configurations

    """

    def __init__(self, element_id, muts, segments, regions_of_interest, vep_reader, signature, config):
        # Input attributes
        self.name = element_id

        self.muts = muts

        self.signature = signature
        self.segments = segments
        self.symbol = self.segments[0].get('SYMBOL', None)

        self.regions_of_interest = regions_of_interest

        # Configuration parameters
        self.sampling_size = config['statistic']['sampling']
        self.sampling_chunk = config['statistic']['sampling_chunk'] * 10**6
#        self.statistic_name = config['statistic']['method']

        # Output attributes
        self.result = {}

    def run(self):
        """
        """

        vep_reader = bgpack_reader.load(bgdata.get_path('bgvep', 'GRCh37', 'most_severe'))

        observed = [(m['POSITION'], m['ALT']) for m in self.muts]

        self.result['muts_count'] = len(self.muts)
        self.result['partitions'] = []
        self.result['in_reg_counts'] = {}
        if len(self.muts) > 0:

            items_to_simulate = []
            items_to_simulate_prob = []

            with vep_reader as reader:

                for region in self.segments:

                    for row in reader.get(region['CHROMOSOME'], region['START'], region['STOP']):

                        pos, conseqs = row

                        ref_triplet = get_ref_triplet(region['CHROMOSOME'], pos - 1)
                        ref = ref_triplet[1]

                        for alt, conseq in zip('ACGT', conseqs):
                            if alt == ref:
                                continue

                            alt_triplet = ref_triplet[0] + alt + ref_triplet[2]

                            items_to_simulate.append((pos, alt))

                            if conseq == 'missense' or conseq == 'missense_variant':
                                if self.signature is None:
                                    items_to_simulate_prob.append(1)
                                else:
                                    items_to_simulate_prob.append(self.signature.get((ref_triplet, alt_triplet), 0.0))
                            else:
                                items_to_simulate_prob.append(0.0)

            if sum(items_to_simulate_prob) == 0:
                logger.warning('Probability of simulation equal to 0 in {}'.format(self.name))
                self.result['partitions'] = []
            else:
                # normalize probs
                items_to_simulate_prob = np.array(items_to_simulate_prob)
                items_to_simulate_prob = items_to_simulate_prob / sum(items_to_simulate_prob)

                muts_count = len(self.muts)
                chunk_count = (self.sampling_size * muts_count) // self.sampling_chunk
                chunk_size = self.sampling_size if chunk_count == 0 else self.sampling_size // chunk_count

                # Calculate sampling parallelization partitions
                self.result['partitions'] = partitions_list(self.sampling_size, chunk_size)
                self.result['sampling'] = {}

                # Run first partition
                first_partition = self.result['partitions'].pop(0)

                indexes=range(len(items_to_simulate))
                background_index = np.random.choice(indexes, size=(first_partition, muts_count), p=items_to_simulate_prob, replace=True)

                # Get the mutations back from indexes

                list_mutations_simulated = []
                for list_muts in background_index:
                    mutations = [items_to_simulate[x] for x in list_muts]
                    list_mutations_simulated.append(mutations)

                # Iterate over the regions of that particular element
                in_reg_counts = {}

                for reg in self.regions_of_interest:
                    name = reg.data
                    start = reg.begin
                    end = reg.end
                    count = 0
                    for mut in observed:
                        if mut[0] >= start and mut[0] <= end:
                            count += 1
                    in_reg_counts[name] = [count]

                for reg in self.regions_of_interest:
                    name = reg.data
                    start = reg.begin
                    end = reg.end
                    for sim in list_mutations_simulated:
                        count = 0
                        for mut in sim:
                            if mut[0] >= start and mut[0] <= end:
                                count += 1
                        in_reg_counts[name].append(count)

                # TODO statitical test and store results

                # Sampling parallelization (if more than one partition)
                if len(self.result['partitions']) > 0:
                    self.result['region_of_interest'] = self.regions_of_interest
                    self.result['muts_count'] = muts_count
                    self.result['simulation_items'] = items_to_simulate
                    self.result['simulation_probs'] = items_to_simulate_prob

            self.result['in_reg_counts'] = in_reg_counts
        self.result['sampling_size'] = self.sampling_size
        self.result['symbol'] = self.symbol
        self.result['observed'] = observed

        return self
