import os
import sys
import csv
import gzip
import signal
import subprocess

from os import path
from .base import Task


class MutsigCvTask(Task):

    KEY = 'mutsigcv'

    MAF_TO_VEP = [
        ("Splice_Site", ["splice_acceptor_variant", "splice_donor_variant", "transcript_ablation", "exon_loss_variant"]),
        ("Nonsense_Mutation", ['stop_gained']),
        ("Nonstop_Mutation", ['stop_lost']),
        ("Translation_Start_Site", ["initiator_codon_variant", "start_lost"]),
        ("Silent", ["incomplete_terminal_codon_variant", "synonymous_variant", "stop_retained_variant", "NMD_transcript_variant"]),
        ("Missense_Mutation", ["missense_variant", "coding_sequence_variant", "conservative_missense_variant", "rare_amino_acid_variant"]),
    ]

    # TODO
    # return "Frame_Shift_Del" if(( $effect eq 'frameshift_variant' or ( $effect eq 'protein_altering_variant' and !$inframe )) and $var_type eq 'DEL' );
    # return "Frame_Shift_Ins" if(( $effect eq 'frameshift_variant' or ( $effect eq 'protein_altering_variant' and !$inframe )) and $var_type eq 'INS' );
    # return "In_Frame_Ins" if( $effect =~ /^(inframe_insertion|disruptive_inframe_insertion)$/ or ( $effect eq 'protein_altering_variant' and $inframe and $var_type eq 'INS' ));
    # return "In_Frame_Del" if( $effect =~ /^(inframe_deletion|disruptive_inframe_deletion)$/ or ( $effect eq 'protein_altering_variant' and $inframe and $var_type eq 'DEL' ));
    # return "Intron" if ( $effect =~ /^(transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)$/ );
    # return "Splice_Region" if( $effect eq 'splice_region_variant' );
    # return "RNA" if( $effect =~ /^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)$/ );
    # return "5'UTR" if( $effect =~ /^(5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant)$/ );
    # return "3'UTR" if( $effect eq '3_prime_UTR_variant' );
    # return "IGR" if( $effect =~ /^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)$/ );
    # return "5'Flank" if( $effect eq 'upstream_gene_variant' );
    # return "3'Flank" if ( $effect eq 'downstream_gene_variant' );

    VEP_TO_MAF_CT = {}
    for m, veps in MAF_TO_VEP:
        for v in veps:
            VEP_TO_MAF_CT[v] = m

    def __init__(self, output_folder, config):

        super().__init__(output_folder, config)

        self.name = None
        self.in_fd = None
        self.in_writer = None
        self.in_file = None
        self.in_skip = False

        self.out_file = None
        self.output_folder = path.join(output_folder, self.KEY)
        os.makedirs(self.output_folder, exist_ok=True)

    def input_start(self):

        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')
            self.in_writer.writerow(["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
                                     "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Consequence", "Variant_Classification"])

    def input_write(self, _, value):

        if not self.in_skip:

            # Hugo_Symbol Chromosome Start_Position End_Position Reference_Allele Tumor_Seq_Allele2 Tumor_Sample_Barcode
            # Variant_Classification Transcript_ID
            if value['CANONICAL'] == 'YES':
                identifier, sample, ref, alt = value['#Uploaded_variation'].split('__')
                chromosome, position = value['Location'].split(':')
                consequence = value['Consequence'].split(',')[0]

                variant_class = self.VEP_TO_MAF_CT.get(consequence, consequence)
                self.in_writer.writerow([
                    value['SYMBOL'],
                    chromosome,
                    position,
                    position,
                    ref,
                    alt,
                    sample,
                    consequence,
                    variant_class
                ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run vep
        if not path.exists(self.out_file):
            cmd = "mkdir -p {2}/tmp && " \
                  "zcat {1} > {2}/tmp/INPUT.maf && " \
                  "singularity exec {0}/mutsigcv.img /opt/MutSigCV_1.4/run_MutSigCV.sh /opt/mcr/v81 " \
                  "{2}/tmp/INPUT.maf {0}/exome_full192.coverage.txt {0}/gene.covariates.txt " \
                  "{2}/tmp/INPUT {0}/mutation_type_dictionary_file.txt {0}/chr_files_hg19 && " \
                  "cp {2}/tmp/INPUT.sig_genes.txt {2}/{3}.out && " \
                  "gzip {2}/{3}.out".format(
                os.environ['MUTSIGCV_DATA'], self.in_file, self.output_folder, self.name)

            try:
                o = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as e:

                # Don't fail if there are not enough mutations just create and empty output
                if "not enough mutations" in e.output.decode():
                    print(e.output.decode())
                    with gzip.open(self.out_file) as fd:
                        fd.write("gene\texpr\treptime\thic\tN_nonsilent\tN_silent\tN_noncoding\tn_nonsilent\tn_silent\tn_noncoding\tnnei\tx\tX\tp\tq\n")
                    sys.exit(0)
                else:
                    print(e.output.decode())
                    sys.exit(e.returncode)
            print(o.decode())

        return self.out_file

    def clean(self):

        if path.exists(self.out_file + ".pid"):
            with open(self.out_file + ".pid", "rt") as fd:
                pid = int(fd.readline().strip())
                os.kill(pid, signal.SIGTERM)
            os.unlink(self.out_file + ".pid")

        if path.exists(self.out_file):
            os.unlink(self.out_file)

        if path.exists(self.in_file):
            os.unlink(self.in_file)
