import os
import csv
import gzip
import itertools
import logging
import numpy as np

from bgreference import refseq
from .base import Filter
from collections import Counter, defaultdict
from intervaltree import IntervalTree
from pyliftover import LiftOver

logger = logging.getLogger(__name__)


class VariantsFilter(Filter):

    KEY = "variants"

    # Minimum cutoff
    MIN_CUTOFF = {"WXS": 1000, "WGS": 10000}
    CHROMOSOMES = set([str(c) for c in range(1, 23)] + ['X', 'Y'])

    # Annotations
    PLATFORM = {}
    GENOMEREF = {}

    
    def __init__(self, parent):
        super().__init__(parent)
        self.liftover = {
            'hg19': LiftOver(os.path.join(os.environ['INTOGEN_DATASETS'], 'preprocess', 'hg19ToHg38.over.chain.gz')),
            'hg38': LiftOver(os.path.join(os.environ['INTOGEN_DATASETS'], 'preprocess', 'hg38ToHg19.over.chain.gz'))
        }

        # Compute chromosomes maximum lenght
        self.chr_maxposition = {'hg38': {}, 'hg19': {}}
        for chr in self.CHROMOSOMES:
            self.chr_maxposition['hg38'][chr] = len(refseq('hg38', chr, 1, -1))
            self.chr_maxposition['hg19'][chr] = len(refseq('hg19', chr, 1, -1))


    @staticmethod
    def __none_to_string(value):
        if value is None:
            return "None"
        return value

    def hypermutators_cutoff(self, snp_per_sample, sequence_type="WXS"):
        vals = list(snp_per_sample.values())
        iqr = np.subtract(*np.percentile(vals, [75, 25]))
        q3 = np.percentile(vals, 75)
        computed_cutoff = (q3 + 1.5 * iqr)
        cutoff = max(self.MIN_CUTOFF[sequence_type], computed_cutoff)
        return cutoff, computed_cutoff, set([k for k, v in snp_per_sample.items() if v > cutoff])

    def sample_stats(self, group_key, group_data):

        donors = {}
        mut_per_sample = {}
        snp_per_sample = {}
        indel_per_sample = {}
        
        for m in self.parent.run(group_key, group_data):

            if group_key not in self.PLATFORM and 'PLATFORM' in m:
                self.PLATFORM[group_key] = m['PLATFORM'].upper()

            if group_key not in self.GENOMEREF and 'GENOMEREF' in m:
                self.GENOMEREF[group_key] = m['GENOMEREF'].lower()

            s = m['SAMPLE']

            if m['ALT_TYPE'] == 'snp':
                snp_per_sample[s] = snp_per_sample.get(s, 0) + 1

            elif m['ALT_TYPE'] == 'indel':
                indel_per_sample[s] = indel_per_sample.get(s, 0) + 1

            mut_per_sample[s] = mut_per_sample.get(s, 0) + 1
            if 'DONOR' in m:
                s = donors.get(m['DONOR'], set())
                s.add(m['SAMPLE'])
                donors[m['DONOR']] = s

        return indel_per_sample, mut_per_sample, snp_per_sample, donors

    @staticmethod
    def generate_kmers(kmer):
        """Create a dictionary with all the possible kmers (trinucleotides or pentanucleotides) and alternates

        Args:
            kmer (int): kmer nucleotides to calculate the signature (3 or 5)

        Returns:
            dict: keys are 'ref kmer --> alt kmer', values are 0
        """
        half_kmer = kmer // 2
        nucleotides = 'ACTG'
        results = set()

        for permutation in itertools.product(nucleotides, repeat=kmer):
            kmer_nucleotide = ''.join(permutation)
            for alt in set(nucleotides).difference(kmer_nucleotide[half_kmer]):
                results.add(
                    (kmer_nucleotide, ''.join([kmer_nucleotide[0:half_kmer], alt, kmer_nucleotide[half_kmer + 1:]]))
                )

        return {'>'.join(key): 0 for key in results}

    def run(self, group_key, group_data):

        group_key = 'None' if group_key is None else group_key

        # To store errors and statistics
        self.stats[group_key] = {}

        # Compute sample and donor statistics
        indel_per_sample, mut_per_sample, snp_per_sample, donors = self.sample_stats(group_key, group_data)
        self.stats[group_key]['samples'] = {
            'count': len(mut_per_sample),
            'mut_per_sample': mut_per_sample,
            'snp_per_sample': snp_per_sample,
            'indel_per_sample': indel_per_sample
        }

        for d in donors:
            donors[d] = list(sorted(donors[d]))

        self.stats[group_key]['donors'] = {self.__none_to_string(d): s for d, s in donors.items()}

        if len(snp_per_sample) < 1:
            self.stats[group_key]['error_no_samples_with_snp'] = "[{}] Any sample has a SNP variant".format(group_key)
            return

        if None in donors:
            self.stats[group_key]['warning_no_donor_id'] = "[{}] There is no DONOR id".format(group_key)

        donors_with_multiple_samples = [d for d, s in donors.items() if len(s) > 1]
        if len(donors_with_multiple_samples) > 0:
            self.stats[group_key]['warning_multiple_samples_per_donor'] = "[{}] {}".format(group_key, donors_with_multiple_samples)

        # We only want to use one sample per donor 
        # this dictionary contains the samples that we'll skip
        multiple_donor_samples = []
        for d in donors_with_multiple_samples:
            multiple_donor_samples += donors[d][1:]
        multiple_donor_samples = set(multiple_donor_samples)        

        sequence_type = self.PLATFORM.get(group_key, "WXS")
        if group_key not in self.PLATFORM:
            self.stats[group_key]['warning_unknown_platform'] = "[{}] There is no PLATFORM annotation, using WXS by default".format(group_key)


        cutoff, theorical_cutoff, hypermutators = self.hypermutators_cutoff(snp_per_sample, sequence_type=sequence_type)
        self.stats[group_key]['hypermutators'] = {
            'cutoff': cutoff,
            'computed_cutoff': theorical_cutoff,
            'hypermutators': list(hypermutators)
        }

        # This cohort genome reference
        genome_ref = self.GENOMEREF.get(group_key, "hg19")
        if group_key not in self.GENOMEREF:
            self.stats[group_key]['warning_unknown_genomeref'] = "[{}] There is no GENOMEREF annotation, using HG19 by default".format(group_key)

        # Load coverage regions tree
        regions_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'preprocess','{}_100bp.coverage.regions.gz'.format(genome_ref))
        coverage_tree = defaultdict(IntervalTree)
        with gzip.open(regions_file, 'rt') as fd:
            reader = csv.reader(fd, delimiter='\t')
            for i, r in enumerate(reader, start=1):
                coverage_tree[r[0]][int(r[1]):(int(r[2]) + 1)] = i

        # Load Somatic Pon file
        somatic_pon_file = os.path.join(
            os.environ['INTOGEN_DATASETS'], 'shared', 'somatic_pon_count_filtered.tsv.gz'
        )
        somatic_pon = set()
        with gzip.open(somatic_pon_file, 'rt') as fd:
            for line in fd:
                line = line.strip().split('\t')
                somatic_pon.add(line)

        # Stats counter
        skip_hypermutators = 0
        skip_multiple_samples_per_donor = 0
        skip_chromosome = 0
        skip_chromosome_names = set()
        skip_coverage = 0
        skip_coverage_positions = []
        skip_same_alt = 0
        skip_n_sequence = 0
        skip_somatic_pon = 0
        skip_mismatch = 0
        skip_duplicated = 0
        skip_no_liftover = 0

        count_before = 0
        count_after = 0
        count_snp = 0
        count_indel = 0
        
        # Read variants
        signature = self.generate_kmers(3)
        variants_by_sample = {}

        for v in self.parent.run(group_key, group_data):
            count_before += 1

            # Force the strand to be positive
            v['STRAND'] = '+'

            if v['REF'] == v['ALT']:
                skip_same_alt += 1
                continue

            # Skip hypermutators and multiple samples per donor           
            if v['SAMPLE'] in hypermutators:
                skip_hypermutators += 1
                continue  

            # Skip multiple samples per donor           
            if v['SAMPLE'] in multiple_donor_samples:
                skip_multiple_samples_per_donor += 1
                continue                

            if v['CHROMOSOME'] not in self.CHROMOSOMES:
                skip_chromosome_names.add(v['CHROMOSOME'])
                skip_chromosome += 1
                continue

            if v['CHROMOSOME'] in coverage_tree:
                if len(coverage_tree[v['CHROMOSOME']][v['POSITION']]) == 0:
                    skip_coverage += 1
                    skip_coverage_positions.append((v['SAMPLE'], v['CHROMOSOME'], v['POSITION']))
                    continue

            # Check duplicates
            var_value = "{}:{}:{}>{}".format(v['CHROMOSOME'], v['POSITION'], v['REF'], v['ALT'])
            if v['SAMPLE'] in variants_by_sample:
                if var_value in variants_by_sample[v['SAMPLE']]:
                    skip_duplicated += 1
                    continue
                else:
                    variants_by_sample[v['SAMPLE']].add(var_value)
            else:
                variants_by_sample[v['SAMPLE']] = {var_value}

            # Skip variants that has an N
            if 'N' in v['ALT'] or 'N' in v['REF']:
                skip_n_sequence += 1
                continue

            # Skip variants that are in the somatic_pon_count_filtered.tsv.gz file
            var_value = (v['CHROMOSOME'], v['POSITION'], v['REF'], v['ALT'])
            if var_value in somatic_pon:
                skip_somatic_pon += 1
                continue

            count_after += 1
            if v['ALT_TYPE'] == 'snp':
                count_snp += 1

                # Compute signature and count mismatch
                ref = refseq(genome_ref, v['CHROMOSOME'], v['POSITION'] - 1, size=3).upper()
                alt = ''.join([ref[0], v['ALT'], ref[2]])
                if ref[1] != v['REF']:
                    skip_mismatch += 1
                    continue
                signature_key = "{}>{}".format(ref, alt)
                signature[signature_key] = signature.get(signature_key, 0) + 1

            elif v['ALT_TYPE'] == 'indel':
                count_indel += 1

            # Add liftover columns
            lo_position = self.liftover[genome_ref].convert_coordinate("chr{}".format(v['CHROMOSOME']), v['POSITION'] - 1, v['STRAND'])
            if lo_position is None or len(lo_position) != 1:
                skip_no_liftover += 1
                continue 

            # Check that the liftover is inside the chromosome
            lo_position = lo_position[0][1] + 1
            if lo_position < 1 or lo_position > self.chr_maxposition[genome_ref][v['CHROMOSOME']]:
                skip_no_liftover += 1
                continue
            
            v['POSITION_HG19'] = v['POSITION'] if genome_ref == 'hg19' else lo_position
            v['POSITION_HG38'] = v['POSITION'] if genome_ref == 'hg38' else lo_position

            yield v

        self.stats[group_key]['signature'] = signature
        signature_count = sum(signature.values())
        if signature_count > 0:
            self.stats[group_key]['probabilities'] = {k: v / signature_count for k, v in signature.items()}
        else:
            self.stats[group_key]['probabilities'] = signature

        self.stats[group_key]['skip'] = {
            'hypermutators': skip_hypermutators,
            'multiple_samples_per_donor': skip_multiple_samples_per_donor,
            'invalid_chromosome': (skip_chromosome, list(skip_chromosome_names)),
            'coverage': (skip_coverage, skip_coverage_positions),
            'same_alt': skip_same_alt,
            'n_sequence': skip_n_sequence,
            'mismatch': skip_mismatch,
            'duplicated': skip_duplicated,
            'noliftover': skip_no_liftover,
            'somatic_pon': skip_somatic_pon
        }

        self.stats[group_key]['count'] = {
            'before': count_before,
            'after': count_after,
            'snp': count_snp,
            'indel': count_indel
        }

        self.stats[group_key]['signature'] = signature

        ratio_mismatch = (skip_mismatch / count_snp) if count_snp > 0 else 0
        if ratio_mismatch > 0.1:
            self.stats[group_key]["error_genome_reference_mismatch"] = "[{}] There are {} of {} genome reference mismatches. More than 10%, skipping this dataset.".format(group_key, skip_mismatch, count_snp)
        elif ratio_mismatch > 0.05:
            self.stats[group_key]["warning_genome_reference_mismatch"] = "[{}] There are {} of {} genome reference mismatches.".format(group_key, skip_mismatch, count_snp)

        if skip_same_alt > 0:
            self.stats[group_key]["warning_same_alternate"] = "[{}] There are {} entries with same reference and alternate".format(group_key, skip_same_alt)

        if count_after == 0:
            self.stats[group_key]["error_no_entries"] = "[{}] There are no variants after filtering".format(group_key)

        if skip_duplicated > 0:
            self.stats[group_key]["warning_duplicated_variants"] = "[{}] There are {} duplicated variants".format(group_key, skip_duplicated)

        if skip_n_sequence > 0:
            self.stats[group_key]["warning_n_sequence"] = "[{}] There are {} variants with a 'N' in the reference or alternate sequence".format(group_key, skip_n_sequence)
