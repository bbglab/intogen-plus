import os
import sys
import csv
import json
import gzip
import subprocess
import logging
import pickle

from os import path, makedirs
from functools import partial
from concurrent.futures import ProcessPoolExecutor


VALID_CONSEQUENCES = {
        "transcript_ablation", "splice_donor_variant", "splice_acceptor_variant", "stop_gained", "frameshift_variant",
        "stop_lost", "initiator_codon_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion",
        "missense_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant",
        "synonymous_variant", "coding_sequence_variant"
    }

logger = logging.getLogger('intogen')

def prepare_task(reader, tasks, item):
    position, (group_key, group_data) = item

    # Initialize task name
    [t.init(group_key) for t in tasks]

    # Input start
    [t.input_start() for t in tasks]

    # Input write
    for i, mut in enumerate(reader.run(group_key, group_data)):
        id = "I{:010d}".format(i)
        for t in tasks:
            t.input_write(id, mut)

    # Input close
    [t.input_end() for t in tasks]

    # Check errors and remove inputs
    errors = [e for e in reader.stats.get(group_key, {}).keys() if e.startswith("error_")]

    if len(errors) > 0:
        for t in tasks:
            if os.path.exists(t.in_file):
                os.unlink(t.in_file)

    return [(t.KEY, group_key) for t in tasks], reader.stats


def prepare_tasks(output, groups, reader, tasks, cores=None):
    func = partial(prepare_task, reader, tasks)

    all_tasks = []
    if cores > 1:
        with ProcessPoolExecutor(cores) as executor:
            for tasks, stats in executor.map(func, enumerate(groups)):
                all_tasks += tasks
                reader.stats.update(stats)
    else:
        for tasks in map(func, enumerate(groups)):
            all_tasks += tasks

    # Store filters stats
    stats_file = os.path.join(output, "filters", "{}.json".format(reader.KEY))
    os.makedirs(os.path.dirname(stats_file), exist_ok=True)

    # Store signatures
    if reader.KEY == 'variants':
        for k, v in reader.stats.items():
            if 'signature' in v and 'probabilities' in v:
                signature_file = os.path.join(output, "signatures", "{}.pickle".format(k))
                os.makedirs(os.path.dirname(signature_file), exist_ok=True)

                with open(signature_file, "wb") as fd:
                    obj = {
                        'counts': {tuple(k_counts.split('>')): v_counts for k_counts, v_counts in v['signature'].items()},
                        'probabilities': {tuple(k_counts.split('>')): v_counts for k_counts, v_counts in v['probabilities'].items()}
                    }
                    pickle.dump(obj, fd)

    while reader.parent is not None:
        try:
            with open(stats_file, "wt") as fd:
                json.dump(reader.stats, fd, indent=4, sort_keys=True)
        except TypeError as e:
            logger.error(e)

        reader = reader.parent

    return all_tasks


def valid_consequence(consequences):
    csq_set = {c for c in consequences.split(',')}
    return len(csq_set & VALID_CONSEQUENCES) > 0

def run_command(cmd):
    try:
        cmd_encoded = cmd.replace("\n", "").replace('\'', '\\\'')
        o = subprocess.check_output(f"bash -c '{cmd_encoded}'", shell=True)
    except subprocess.CalledProcessError as e:
        print(e.output.decode())
        sys.exit(e.returncode)

    print(o.decode())


class Task:

    KEY = None
    INPUT_HEADER = None

    def __init__(self, output_folder, output_project=None, output_work=None):

        self.name = None
        self.in_fd = None
        self.in_file = None
        self.process = None
        self.out_file = None
        self.signatures_file = None
        self.output_project = output_project
        self.output_work = output_work

        self.output_folder = path.join(output_folder, self.KEY)
        makedirs(self.output_folder, exist_ok=True)

        # Default command line and datasets folder
        self.cmdline = os.path.join(os.environ['INTOGEN_CONTAINERS'], f"{self.KEY}.simg")
        self.datasets = os.path.join(os.environ['INTOGEN_DATASETS'], f"{self.KEY}")

        # Load selected transcripts
        with open(os.path.join(os.environ['INTOGEN_DATASETS'], 'shared', 'ensembl_canonical_transcripts.tsv')) as fd:
            self.transcripts = set([r[1] for r in csv.reader(fd, delimiter='\t')])

    def input_start(self):
        self.in_fd = gzip.open(self.in_file, 'wt')
        self.in_writer = csv.writer(self.in_fd, delimiter='\t')
        if self.INPUT_HEADER:
            self.in_writer.writerow(self.INPUT_HEADER)

    def input_write(self, _, value):
        pass

    def write_row(self, row):
        self.in_writer.writerow(row)

    def input_end(self):
        if self.in_fd:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):
        raise NotImplementedError("Class must implement method run()")

    def init(self, name):

        if name is None:
            name = "None"

        if path.isabs(name):
            self.in_file = name
            self.name = path.basename(name).split('.')[0]
        else:
            self.name = name
            self.in_file = path.join(self.output_folder, "{}.in.gz".format(name))

        self.process = None
        self.out_file = path.join(self.output_folder, "{}.out.gz".format(self.name))

        if self.output_project is not None:
            self.signatures_file = os.path.join(self.output_project, 'signatures', f'{self.name}.json')

    def __repr__(self):
        return "{} '{}'".format(self.KEY, self.name)
