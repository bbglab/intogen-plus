import csv
import gzip
import os

from os.path import join, exists
from oncodriveclust.analysis import OncodriveClustAnalysis, get_cluster_coordinates_output, sort_matrix, add_fdr

# Consequence type
SYNONYMOUS = {"synonymous_variant"}
NON_SYNONYMOUS = {"missense_variant", "stop_gained", "stop_lost", "initiator_codon_variant",
                  "incomplete_terminal_codon_variant", "splice_donor_variant", "splice_acceptor_variant"}


class OncodriveClustTask:

    def __init__(self, output_folder):
        self.name = None

        self.min_gene_mutations = 5
        self.intraclust_max_distance = 5
        self.binom_noise_p = 0.01

        self.cds_len = None
        self.non_syn = None
        self.syn = None

        self.out_file = None
        self.output_folder = join(output_folder, "oncodriveclust")
        os.makedirs(self.output_folder, exist_ok=True)

    def __repr__(self):
        return "OncodriveClust '{}'".format(self.name)

    def init(self, name):
        self.name = name
        self.out_file = join(self.output_folder, "{}.out.gz".format(name))

    def input_start(self):
        self.non_syn = {}
        self.syn = {}
        self.cds_len = {}

    @staticmethod
    def _update_dict(data, gene, position):

        if gene not in data:
            data[gene] = {}

        if position not in data[gene]:
            data[gene][position] = 0

        data[gene][position] += 1

    def input_write(self, _, value):
        if value['CANONICAL'] == 'YES':
            try:
                gene = value['Gene']
                position = int(value['Protein_position'])
                consequences = {c.strip() for c in value['Consequence'].split(",")}

                if consequences.intersection(SYNONYMOUS):
                    self._update_dict(self.syn, gene, position)
                    self.cds_len[gene] = 1234

                elif consequences.intersection(NON_SYNONYMOUS):
                    self._update_dict(self.non_syn, gene, position)
                    self.cds_len[gene] = 1234
            except ValueError:
                pass

    def input_end(self):
        pass

    def run(self):

        analysis = OncodriveClustAnalysis()

        non_syn_accum_mut_pos, non_syn_cluster_coordinates, syn_gene_cluster_scores,\
        non_syn_gene_cluster_scores, non_syn_cluster_muts, non_syn_gene_cluster_scores_external_z = \
            analysis.run(
            self.non_syn, self.syn,
            self.min_gene_mutations, self.cds_len, self.intraclust_max_distance, self.binom_noise_p
        )

        header = ['GENE', 'GENE_LEN', 'GENE_NUM_MUTS', 'MUTS_IN_CLUST', 'NUM_CLUSTERS', 'CLUST_COORDS', 'GENE_SCORE', 'ZSCORE', 'PVALUE', 'QVALUE']

        with gzip.open(self.out_file, "wt") as outf:
            writer = csv.writer(outf, delimiter='\t')

            writer.writerow(header)

            output = []
            for gene in non_syn_gene_cluster_scores:
                gene_len = int(self.cds_len[gene]) / 3
                gene_muts = sum([non_syn_accum_mut_pos[gene][pos] for pos in list(non_syn_accum_mut_pos[gene].keys())])
                num_clusters = len(list(non_syn_cluster_coordinates[gene].keys()))

                if num_clusters > 0:
                    muts_clusters = sum(non_syn_cluster_muts[gene].values())
                    gene_additive_cluster_score = non_syn_gene_cluster_scores[gene]
                    zscore = non_syn_gene_cluster_scores_external_z[
                        gene] if gene in non_syn_gene_cluster_scores_external_z else 'NA'
                    row = [gene, gene_len, gene_muts, muts_clusters, num_clusters]
                    cluster_coordinates = get_cluster_coordinates_output(gene, non_syn_cluster_coordinates,
                                                                         non_syn_cluster_muts)
                    row += [cluster_coordinates]
                    row += [gene_additive_cluster_score, zscore]
                    output.append(row)

            # Sort by p-value
            output = sort_matrix(output, -1)

            # Compute q-value
            output = add_fdr(output, -1)

            for row in output:
                writer.writerow(row)

    def clean(self):

        if exists(self.out_file):
            os.unlink(self.out_file)
